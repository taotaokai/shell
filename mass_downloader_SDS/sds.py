# -*- coding: utf-8 -*-
"""
obspy.clients.filesystem.sds - read support for SeisComP Data Structure
=======================================================================
This module provides read support for data stored locally in a SeisComP Data
Structure (SDS) directory structure.

The directory and file layout of SDS is defined as::

    <SDSdir>/YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY

These fields are defined by SDS as follows::

    SDSdir :  arbitrary base directory
    YEAR   :  4 digit year
    NET    :  Network code/identifier, up to 8 characters, no spaces
    STA    :  Station code/identifier, up to 8 characters, no spaces
    CHAN   :  Channel code/identifier, up to 8 characters, no spaces
    TYPE   :  1 characters indicating the data type, recommended types are:
               'D' - Waveform data
               'E' - Detection data
               'L' - Log data
               'T' - Timing data
               'C' - Calibration data
               'R' - Response data
               'O' - Opaque data
    LOC    :  Location identifier, up to 8 characters, no spaces
    DAY    :  3 digit day of year, padded with zeros

See https://www.seiscomp3.org/wiki/doc/applications/slarchive/SDS.

:copyright:
    The ObsPy Development Team (devs@obspy.org)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""
import glob
import os
import re
import warnings
import logging
from datetime import timedelta

import numpy as np

from obspy import Stream, read, UTCDateTime
from obspy.core.stream import _headonly_warning_msg
from obspy.core.util.misc import BAND_CODE
from obspy.io.mseed import ObsPyMSEEDFilesizeTooSmallError
from obspy.io.mseed.util import get_record_information

SDS_FMTSTR = os.path.join(
    "{year}", "{network}", "{station}", "{channel}.{sds_type}",
    "{network}.{station}.{location}.{channel}.{sds_type}.{year}.{doy:03d}")
FORMAT_STR_PLACEHOLDER_REGEX = r"{(\w+?)?([!:].*?)?}"  # this also works: r'{(\w+):?\w*}'


class Client(object):
    """
    Request client for SeisComP Data Structure archive on local filesystem.

    For details see the :meth:`~obspy.clients.filesystem.sds.Client.__init__()`
    method.
    """
    FMTSTR = SDS_FMTSTR

    def __init__(self, sds_root, sds_type="D", format="MSEED",
                 fileborder_seconds=30, fileborder_samples=5000, logger=None):
        """
        Initialize a SDS local filesystem client.

        >>> from obspy.clients.filesystem.sds import Client
        >>> client = Client("/my/SDS/archive/root")  # doctest: +SKIP

        :type sds_root: str
        :param sds_root: Root directory of SDS archive.
        :type sds_type: str
        :param sds_type: SDS data type identifier, one single character. Types
            recommended by the SDS definition are: 'D' - Waveform data,
            'E' - Detection data, 'L' - Log data, 'T' - Timing data,
            'C' - Calibration data, 'R' - Response data, 'O' - Opaque data. Can
            also be wildcarded by setting to ``?`` or ``*``.
        :type format: str
        :param format: File format the data is stored in, see
            :func:`~obspy.core.stream.read()` for a list of file formats
            supported by ObsPy. Usually, SDS archives are stored in "MSEED"
            format. Can be set to ``None`` for file format autodetection
            (slowing down the reading).
        :type fileborder_seconds: float
        :param fileborder_seconds: Defines in which cases the client checks the
            previous/next daily file for the requested data (data in SDS
            archives usually spill over on the day break, at least for a few
            seconds). For example setting ``fileborder_seconds=30`` means that
            for a request with ``starttime`` at ``00:00:29`` UTC (or an
            ``endtime`` at ``23:59:31``), the previous daily file is also
            checked, if it contains matching data. The maximum of both
            ``fileborder_seconds`` and ``fileborder_samples`` is used when
            determining if previous/next day should be checked for data.
        :type fileborder_samples: int
        :param fileborder_samples: Similar to ``fileborder_seconds``. The given
            number of samples is converted to seconds by mapping the band
            code of the requested channel to sampling frequency. The maximum of
            both ``fileborder_seconds`` and ``fileborder_samples`` is used when
            determining if previous/next day should be checked for data.
        """
        if not os.path.isdir(sds_root):
            msg = ("SDS root is not a local directory: " + sds_root)
            raise IOError(msg)
        self.sds_root = sds_root
        self.sds_type = sds_type
        self.format = format and format.upper()
        self.fileborder_seconds = fileborder_seconds
        self.fileborder_samples = fileborder_samples
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger(__name__)

    def get_waveforms(self, network, station, location, channel, starttime,
                      endtime, merge=-1, sds_type=None, **kwargs):
        """
        Read data from a local SeisComP Data Structure (SDS) directory tree.

        >>> from obspy import UTCDateTime
        >>> t = UTCDateTime("2015-10-12T12")
        >>> st = client.get_waveforms("IU", "ANMO", "*", "HH?", t, t+30)
        ... # doctest: +SKIP

        :type network: str
        :param network: Network code of requested data (e.g. "IU").
            Wildcards '*' and '?' are supported.
        :type station: str
        :param station: Station code of requested data (e.g. "ANMO").
            Wildcards '*' and '?' are supported.
        :type location: str
        :param location: Location code of requested data (e.g. "").
            Wildcards '*' and '?' are supported.
        :type channel: str
        :param channel: Channel code of requested data (e.g. "HHZ").
            Wildcards '*' and '?' are supported.
        :type starttime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param starttime: Start of requested time window.
        :type endtime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param endtime: End of requested time window.
        :type merge: int or None
        :param merge: Specifies, which merge operation should be performed
            on the stream before returning the data. Default (``-1``) means
            only a conservative cleanup merge is performed to merge seamless
            traces (e.g. when reading across day boundaries). See
            :meth:`Stream.merge(...) <obspy.core.stream.Stream.merge>` for
            details. If set to ``None`` (or ``False``) no merge operation at
            all will be performed.
        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :param kwargs: Additional kwargs that get passed on to
            :func:`~obspy.core.stream.read` internally, mostly for internal
            low-level purposes used by other methods.
        :rtype: :class:`~obspy.core.stream.Stream`
        """
        if starttime >= endtime:
            msg = ("'endtime' must be after 'starttime'.")
            raise ValueError(msg)
        sds_type = sds_type or self.sds_type

        seed_pattern = ".".join((network, station, location, channel))

        st = Stream()
        full_paths = self._get_filenames(
            network=network, station=station, location=location,
            channel=channel, starttime=starttime, endtime=endtime,
            sds_type=sds_type, check_exist=True)
        for full_path in full_paths:
            try:
                # for headonly read, no trim/merge operation is allowed
                if kwargs.get("headonly", False):
                    st += read(full_path, format=self.format, sourcename=seed_pattern, **kwargs)
                else:
                    st += read(full_path, format=self.format, starttime=starttime,
                                endtime=endtime, sourcename=seed_pattern, **kwargs)
            except ObsPyMSEEDFilesizeTooSmallError:
                # just ignore small MSEED files, in use cases working with
                # near-realtime data these are usually just being created right
                # at request time, e.g. when fetching current data right after
                # midnight
                continue
            except Exception as e:
                msg = f'failed to read mseed file in SDS archive: {full_path} (Error: {e})'
                self.logger.warning(msg)
                continue

        # make sure we only have the desired data, just in case the file
        # contents do not match the expected SEED id
        st = st.select(network=network, station=station, location=location,
                       channel=channel)

        # avoid trim/merge operations when we do a headonly read for
        # `_get_availability_percentage()`
        if kwargs.get("_no_trim_or_merge", False) or kwargs.get("headonly", False):
            return st

        st.trim(starttime, endtime)
        if merge is None or merge is False:
            pass
        else:
            st.merge(merge)
        return st

    def get_waveforms_bulk(self, bulk):
        """
        Reads bulk data from a local SeisComP Data Structure (SDS) directory
        tree.

        Returns a stream object with with the data requested from the local
        directory.

        :type bulk: list[tuple]
        :param bulk: Information about the requested data.
        """
        st = Stream()
        for bulk_string in bulk:
            st += self.get_waveforms(*bulk_string)
        return st

    def _get_filenames(self, network, station, location, channel, starttime,
                       endtime, sds_type=None, check_exist=False):
        """
        Get list of filenames for certain waveform and time span.

        :type network: str
        :param network: Network code of requested data (e.g. "IU").
        :type station: str
        :param station: Station code of requested data (e.g. "ANMO").
        :type location: str
        :param location: Location code of requested data (e.g. "").
        :type channel: str
        :param channel: Channel code of requested data (e.g. "HHZ").
        :type starttime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param starttime: Start of time span.
        :type endtime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param endtime: End of time span.
        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :type check_exist: bool
        :param check_exist: Check if file exists.
        :rtype: list[str]
        """
        sds_type = sds_type or self.sds_type
        # SDS has data sometimes in adjacent days, so also try to read the
        # requested data from those files. Usually this is only a few seconds
        # of data after midnight, but for now we play safe here to catch all
        # requested data (and with MiniSEED - the usual SDS file format - we
        # can use starttime/endtime kwargs anyway to read only desired parts).
        year_doy = set()
        # determine how far before starttime/after endtime we should check
        # other dayfiles for the data
        t_buffer = self.fileborder_samples / BAND_CODE.get(channel[:1], 20.0)
        t_buffer = max(t_buffer, self.fileborder_seconds)
        t = starttime - t_buffer
        t_max = endtime + t_buffer
        # make a list of year/doy combinations that covers the whole requested
        # time window (plus day before and day after)
        while t < t_max:
            year_doy.add((t.year, t.julday))
            t += timedelta(days=1)
        year_doy.add((t_max.year, t_max.julday))

        # full_paths = set()
        full_paths = []
        for year, doy in year_doy:
            filename = self.FMTSTR.format(
                network=network, station=station, location=location,
                channel=channel, year=year, doy=doy, sds_type=sds_type)
            full_path = os.path.join(self.sds_root, filename)
            if (not check_exist) or (check_exist and os.path.isfile(full_path)):
                full_paths.append(os.path.join(self.sds_root, filename))

        return full_paths

    def _get_filename(self, network, station, location, channel, time,
                      sds_type=None):
        """
        Get filename for certain waveform.

        :type network: str
        :param network: Network code of requested data (e.g. "IU").
        :type station: str
        :param station: Station code of requested data (e.g. "ANMO").
        :type location: str
        :param location: Location code of requested data (e.g. "").
        :type channel: str
        :param channel: Channel code of requested data (e.g. "HHZ").
        :type time: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param time: Time of interest.
        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :rtype: str
        """
        sds_type = sds_type or self.sds_type
        filename = self.FMTSTR.format(
            network=network, station=station, location=location,
            channel=channel, year=time.year, doy=time.julday,
            sds_type=sds_type)
        return os.path.join(self.sds_root, filename)

    def get_availability_percentage(self, network, station, location, channel,
                                    starttime, endtime, sds_type=None):
        """
        Get percentage of available data.

        :type network: str
        :param network: Network code of requested data (e.g. "IU").
        :type station: str
        :param station: Station code of requested data (e.g. "ANMO").
        :type location: str
        :param location: Location code of requested data (e.g. "").
        :type channel: str
        :param channel: Channel code of requested data (e.g. "HHZ").
        :type starttime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param starttime: Start of requested time window.
        :type endtime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param endtime: End of requested time window.
        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :rtype: tuple(float, int)
        :returns: 3-tuple of percentage of available data (``0.0`` to ``1.0``)
            , number of gaps and list of gaps.
        """
        if starttime >= endtime:
            msg = ("'endtime' must be after 'starttime'.")
            raise ValueError(msg)
        sds_type = sds_type or self.sds_type
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", _headonly_warning_msg, UserWarning,
                "obspy.core.stream")
            st = self.get_waveforms(network, station, location, channel,
                                    starttime, endtime, sds_type=sds_type,
                                    headonly=True, _no_trim_or_merge=True)
        # even if the warning was silently caught and not shown it gets
        # registered in the __warningregistry__ and will not be shown
        # subsequently in a place were it's not caught
        # see https://bugs.python.org/issue4180
        # see e.g. http://blog.ionelmc.ro/2013/06/26/testing-python-warnings/
        try:
            from obspy.core.stream import __warningregistry__ as \
                stream_warningregistry
        except ImportError:
            # import error means no warning has been issued from
            # obspy.core.stream before, so nothing to do.
            pass
        else:
            for key in list(stream_warningregistry.keys()):
                if key[0] == _headonly_warning_msg:
                    stream_warningregistry.pop(key)
        # st.sort(keys=['starttime', 'endtime'])
        # st.traces = [tr for tr in st
        #              if not (tr.stats.endtime < starttime or
        #                      tr.stats.starttime > endtime)]

        # if not st:
        #     return (0, 1)

        # # sum up gaps in the middle
        # gaps = [gap[6] for gap in st.get_gaps()]
        # gap_sum = np.sum(gaps)
        # gap_count = len(gaps)
        # # check if we have a gap at start or end
        # earliest = min([tr.stats.starttime for tr in st])
        # latest = max([tr.stats.endtime for tr in st])
        # if earliest > starttime:
        #     gap_sum += earliest - starttime
        #     gap_count += 1
        # if latest < endtime:
        #     gap_sum += endtime - latest
        #     gap_count += 1

        # return (1 - (gap_sum / total_duration), gap_count)

        gaps = []
        if st:
            existing_intervals = [(tr.stats.starttime, tr.stats.endtime) for tr in st]
            existing_intervals = sorted(existing_intervals, key=lambda x: x[0])
            dt = min([tr.stats.delta for tr in st])
            b = starttime
            for w in existing_intervals:
                if w[1] < starttime: continue
                if w[0] > endtime: break
                if w[0] > b+dt:
                    gaps.append((b, min(w[0], endtime)))
                b = w[1]
                if b > endtime: break
            if b < endtime:
                gaps.append((b, endtime))
        else:
            gaps.append((starttime, endtime))

        gap_count = len(gaps)
        gap_sum = np.sum([w[1] - w[0] for w in gaps])
        total_duration = endtime - starttime

        return (1 - (gap_sum / total_duration), gap_count, gaps)

    def _get_current_endtime(self, network, station, location, channel,
                             sds_type=None, stop_time=None,
                             check_has_no_data=True):
        """
        Get time of last sample for given stream.

        ``None`` is returned if no data at all is encountered when going
        backwards until `stop_time` (defaults to Jan 1st 1950).

        :type network: str
        :param network: Network code of requested data (e.g. "IU").
        :type station: str
        :param station: Station code of requested data (e.g. "ANMO").
        :type location: str
        :param location: Location code of requested data (e.g. "").
        :type channel: str
        :param channel: Channel code of requested data (e.g. "HHZ").
        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :type stop_time: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param stop_time: Time at which the search for data is stopped and
            ``None`` is returned. If not specified, stops at ``1950-01-01T00``.
        :type check_has_no_data: bool
        :param check_has_no_data: Whether to perform a check with
            :meth:`has_data` if any data is available at all for the given SEED
            ID. Turns out that this check can take a long time on slow
            filesystems (e.g. NFS), so it might actually make the latency check
            take longer than necessary, so deactivating it might be worth a try
            if speed issues occur.
        :rtype: :class:`~obspy.core.utcdatetime.UTCDateTime` or ``None``
        """
        sds_type = sds_type or self.sds_type

        seed_pattern = ".".join((network, station, location, channel))

        if check_has_no_data and not self.has_data(
                network=network, station=station, location=location,
                channel=channel, sds_type=sds_type):
            return None

        stop_time = stop_time or UTCDateTime(1950, 1, 1)
        st = None
        time = UTCDateTime()

        while not st:
            if time < stop_time:
                return None
            filename = self._get_filename(
                network=network, station=station, location=location,
                channel=channel, time=time, sds_type=sds_type)
            if os.path.isfile(filename):
                try:
                    st = read(filename, format=self.format, headonly=True,
                              sourcename=seed_pattern)
                except ObsPyMSEEDFilesizeTooSmallError:
                    # just ignore small MSEED files, in use cases working with
                    # near-realtime data these are usually just being created
                    # right at request time, e.g. when fetching current data
                    # right after midnight
                    st = None
                else:
                    st = st.select(network=network, station=station,
                                   location=location, channel=channel)
                if st:
                    break
            time -= 24 * 3600

        return max([tr.stats.endtime for tr in st])

    def get_latency(self, network, station, location, channel,
                    sds_type=None, stop_time=None, check_has_no_data=True):
        """
        Get latency for given stream, i.e. difference of current time and
        latest available data for stream in SDS archive. ``None`` is returned
        if no data at all is encountered when going backwards until
        `stop_time` (defaults to Jan 1st 1950).

        :type network: str
        :param network: Network code of requested data (e.g. "IU").
        :type station: str
        :param station: Station code of requested data (e.g. "ANMO").
        :type location: str
        :param location: Location code of requested data (e.g. "").
        :type channel: str
        :param channel: Channel code of requested data (e.g. "HHZ").
        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :type stop_time: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param stop_time: Time at which the search for data is stopped and
            ``None`` is returned. If not specified, stops at ``1950-01-01T00``.
        :type check_has_no_data: bool
        :param check_has_no_data: Whether to perform a check with
            :meth:`has_data` if any data is available at all for the given SEED
            ID. Turns out that this check can take a long time on slow
            filesystems (e.g. NFS), so it might actually make the latency check
            take longer than necessary, so deactivating it might be worth a try
            if speed issues occur.
        :rtype: float or ``None``
        :returns: Latency in seconds or ``None`` if no data was encountered
            from current time backwards until ``stop_time``.
        """
        endtime = self._get_current_endtime(
            network=network, station=station, location=location,
            channel=channel, sds_type=sds_type, stop_time=stop_time,
            check_has_no_data=check_has_no_data)

        if endtime is None:
            return endtime

        return UTCDateTime() - endtime

    def has_data(self, network, station, location, channel, sds_type=None):
        """
        Check if specified stream has any data.

        Actually just checks whether a file is encountered in a folder that is
        expected to contain data.

        :type network: str
        :param network: Network code of requested data (e.g. "IU").
        :type station: str
        :param station: Station code of requested data (e.g. "ANMO").
        :type location: str
        :param location: Location code of requested data (e.g. "").
        :type channel: str
        :param channel: Channel code of requested data (e.g. "HHZ").
        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :rtype: bool
        """
        sds_type = sds_type or self.sds_type

        pattern = re.sub(
            FORMAT_STR_PLACEHOLDER_REGEX,
            _wildcarded_except(["network", "station", "location", "channel",
                                "sds_type"]),
            self.FMTSTR)
        pattern = pattern.format(
            network=network, station=station, location=location,
            channel=channel, sds_type=sds_type)
        pattern = os.path.join(self.sds_root, pattern)
        if glob.glob(pattern):
            return True
        else:
            return False

    def get_all_nslc(self, sds_type=None, datetime=None):
        """
        Return information on what streams are included in archive.

        Note that this can be very slow on network file systems because every
        single file has to be touched (because available location codes can not
        be discovered from folder structure alone).

        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :type datetime: :class:`~obspy.core.utcdatetime.UTCDateTime`
        :param datetime: Only return all streams that have data at specified
            time (checks if file exists that should have the data, i.e. streams
            might be returned that have data on the same day but not at exactly
            this point in time).
        :rtype: list
        :returns: List of (network, station, location, channel) 4-tuples of all
            available streams in archive.
        """
        sds_type = sds_type or self.sds_type
        result = set()
        # wildcarded pattern to match all files of interest
        if datetime is None:
            pattern = re.sub(
                FORMAT_STR_PLACEHOLDER_REGEX,
                _wildcarded_except(["sds_type"]),
                self.FMTSTR).format(sds_type=sds_type)
            pattern = os.path.join(self.sds_root, pattern)
        else:
            pattern = self._get_filename("*", "*", "*", "*", datetime)
        all_files = glob.glob(pattern)
        # set up inverse regex to extract kwargs/values from full paths
        pattern_ = os.path.join(self.sds_root, self.FMTSTR)
        group_map = {i: groups[0] for i, groups in
                     enumerate(re.findall(FORMAT_STR_PLACEHOLDER_REGEX,
                                          pattern_))}
        for file_ in all_files:
            dict_ = _parse_path_to_dict(file_, pattern_, group_map)
            try:
                network = dict_["network"]
                station = dict_["station"]
                location = dict_["location"]
                channel = dict_["channel"]
            except Exception as e:
                msg = (
                    "Failed to extract key from pattern '{}' in path "
                    "'{}': {}").format(pattern_, file_, e)
                self.logger.warning(msg)
                continue
            result.add((network, station, location, channel))
        return sorted(result)

    def get_all_stations(self, sds_type=None):
        """
        Return information on what stations are included in archive.

        This method assumes that network/station combinations can be discovered
        from the folder structure alone (as opposed to the filenames).

        :type sds_type: str
        :param sds_type: Override SDS data type identifier that was specified
            during client initialization.
        :rtype: list
        :returns: List of (network, station) 2-tuples of all available stations
            in archive.
        """
        sds_type = sds_type or self.sds_type
        result = set()
        # wildcarded pattern to match all files of interest
        fmtstr = os.path.dirname(self.FMTSTR)
        pattern = re.sub(
            FORMAT_STR_PLACEHOLDER_REGEX,
            _wildcarded_except(["sds_type"]),
            fmtstr).format(sds_type=sds_type)
        pattern = os.path.join(self.sds_root, pattern)
        all_files = glob.glob(pattern)
        # set up inverse regex to extract kwargs/values from full paths
        pattern_ = os.path.join(self.sds_root, fmtstr)
        group_map = {i: groups[0] for i, groups in
                     enumerate(re.findall(FORMAT_STR_PLACEHOLDER_REGEX,
                                          pattern_))}
        for file_ in all_files:
            dict_ = _parse_path_to_dict(file_, pattern_, group_map)
            try:
                network = dict_["network"]
                station = dict_["station"]
            except (TypeError, KeyError) as e:
                if isinstance(e, TypeError) and dict_ is not None:
                    raise
                msg = (
                    "Failed to extract key from pattern '{}' in path "
                    "'{}': {}").format(pattern_, file_, e)
                self.logger.warning(msg)
                continue
            result.add((network, station))
        return sorted(result)

    def insert_mseed_file(self, mseed_file, sds_type=None, timeout=10):
        from filelock import Timeout, FileLock

        def _parse_mseed(mseed_file):
            mseed_records = []
            file_size = os.path.getsize(mseed_file)
            with open(mseed_file, "rb") as fh:
                try:
                    while True:
                        if fh.tell() >= (file_size - 256):
                            break
                        info = get_record_information(fh)
                        length = info["record_length"]
                        rec = fh.read(length)
                        mseed_records.append((info, rec))
                except Exception as err:
                    msg = f'failed to read {mseed_file}. (Error: {err})'
                    self.logger.warning(msg)
            return mseed_records

        # group mseed records by the corresponding SDS file path
        # groups[db_path] = [(starttime, samp_rate, npts, msr_bytes), ...]
        def _group_msr_for_SDS(mseed_file, sds_type):
            msr_groups = {}
            mseed_records = _parse_mseed(mseed_file)
            for info, rec in mseed_records:
                # strip().strip('\0'): remove whitespaces and null's
                # otherwise, invalid file path containing null might occur 
                net = info['network'].strip().strip('\0')
                sta = info['station'].strip().strip('\0')
                loc = info['location'].strip().strip('\0')
                cha = info['channel'].strip().strip('\0')
                starttime = info['starttime']
                db_path = self._get_filename(net, sta, loc, cha, starttime, sds_type)
                # only keep records with integer sampling rate? 
                samp_rate = info['samp_rate']
                # integer_sr = round(samp_rate)
                # decimal_sr = abs(samp_rate - integer_sr)
                # if decimal_sr != 0:
                #     msg = f'non-integer sampling rate: {db_path} {samp_rate}'
                #     self.logger.warning(msg)
                #     # continue
                npts = info['npts']
                if db_path not in msr_groups:
                    msr_groups[db_path] = [(starttime, samp_rate, npts, rec)]
                else:
                    msr_groups[db_path].append((starttime, samp_rate, npts, rec))
            return msr_groups

        if not os.path.exists(mseed_file):
            msg = f'{mseed_file} does not exit.'
            self.logger.warning(msg)
            return
        elif not os.path.isfile(mseed_file):
            msg = f'{mseed_file} is not a file.'
            self.logger.warning(msg)
            return
        try:
            st = read(mseed_file, format='MSEED', headonly=True)
        except Exception as err:
            msg = f'{mseed_file} is not a mseed file. Error: {err})'
            self.logger.warning(msg)
            return

        sds_type = sds_type or self.sds_type
        msr_groups = _group_msr_for_SDS(mseed_file, sds_type=sds_type)
        for db_path, mseed_records in msr_groups.items():
            try:
                os.makedirs(os.path.dirname(db_path), exist_ok=True)
            except Exception as e:
                msg = f'os.makedirs failed on {db_path.encode()}'
                self.logger.warning(msg)
                continue

            dirname = os.path.dirname(db_path)
            filename = os.path.basename(db_path)
            lock_path = os.path.join(dirname, f'.{filename}.lock')
            lock = FileLock(lock_path, timeout=timeout)
            try:
                with lock:
                    db_records = []
                    if os.path.exists(db_path):
                        try:
                            db_msr_groups = _group_msr_for_SDS(db_path, sds_type=sds_type)
                            if db_path in db_msr_groups:
                                db_records = db_msr_groups[db_path] # should have only one group of mseed records
                        except Exception as e:
                            msg = f'failed to read {db_path} (ERROR: {e})'
                            self.logger.warning(msg)
                            # os.remove(db_path)
                    db_records.extend(mseed_records)

                    # remove duplicated mseed record and update db file
                    updated_records = {}
                    for starttime, samp_rate, npts, rec in sorted(db_records): #, key=lambda x: (x[0], x[1], x[2])):
                        k = (starttime.ns, samp_rate, npts)
                        if k not in updated_records:
                            updated_records[k] = rec # mseed record bytes
                        else:
                            if rec != updated_records[k]:
                                msg = (f'overlap mseed record with different content between ' 
                                       f'{db_path} and {mseed_file} at {starttime}')
                                self.logger.warning(msg)
                            # else:
                            #     warnings.warn(f'duplicated mseed record in {db_path} at {starttime}')

                    mseed_bytes = b''
                    for k, v in sorted(updated_records.items()):
                        mseed_bytes += v
                    with open(db_path, 'wb') as f:
                        f.write(mseed_bytes)
            except Timeout:
                msg = f'cannot insert {mseed_file} since another process currently holds the lock of {db_path}'
                self.logger.warning(msg)


def _wildcarded_except(exclude=[]):
    """
    Function factory for :mod:`re` ``repl`` functions used in :func:`re.sub`,
    replacing all format string place holders with ``*`` wildcards, except
    named fields as specified in ``exclude``.
    """
    def _wildcarded(match):
        if match.group(1) in exclude:
            return match.group(0)
        return "*"
    return _wildcarded


def _parse_path_to_dict(path, pattern, group_map):
    # escape special regex characters "." and "\"
    # in principle we should escape all special characters in Python regex:
    #    . ^ $ * + ? { } [ ] \ | ( )
    # and also only replace them if outside of format string placeholders..
    regex = re.sub(r'([\.\\])', r'\\\1', pattern)
    # replace each format string placeholder with a regex group, matching
    # alphanumerics. append end-of-line otherwise the last non-greedy match
    # doesn't catch anything if it's at the end of the regex
    regex = re.sub(FORMAT_STR_PLACEHOLDER_REGEX, r'(\\w*?)', regex) + "$"
    match = re.match(regex, path)
    if match is None:
        return None
    result = {}
    for i, value in enumerate(match.groups()):
        key = group_map[i]
        if key in result:
            if result[key] != value:
                msg = (
                    "Failed to parse path '{}': Mismatching information "
                    "for key '{}' from pattern '{}'.").format(
                        path, key, pattern)
                warnings.warn(msg)
                return None
        result[key] = value
    return result


if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)
