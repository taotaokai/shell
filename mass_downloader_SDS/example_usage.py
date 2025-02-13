import os, sys
from obspy import UTCDateTime
from obspy.clients.fdsn.mass_downloader import CircularDomain, RectangularDomain, Restrictions
from mass_downloader_SDS import MassDownloader

def download_data(domain, restrictions, mseed_dir, stationxml_dir, log_dir, max_retry_num=5, providers=None):
    ERROR_LIST = [
            "ERROR",
            # "No data available for request",
            # "HTTP Status code: 204",
            # "timed out",
            ]
    # max_retry_num = 5
    current_retry_num = 0
    # Re-launch downloading if error occurs
    while True:
        print(f'\n# retry num: {current_retry_num}\n')

        timestamp = UTCDateTime.now().strftime('%Y%m%dT%H%M%S')
        log_file = f'{log_dir}/mass_downloader_{timestamp}.log'
        with open(log_file, 'w') as f:
            f.write(f'\n\n###### [{timestamp}] retry num: {current_retry_num}\n\n')

        # No specified providers will result in all known ones being queried.
        mdl = MassDownloader(providers=providers, log_file=log_file)
        mdl.download(domain, restrictions, mseed_storage=mseed_dir, stationxml_storage=stationxml_dir)

        retry = False
        with open(log_file, 'r') as f:
            lines = f.readlines()
        for l in lines:
            for err in ERROR_LIST:
                if err in l:
                    retry = True
                    break
            if retry: break
        if not retry: break

        current_retry_num = current_retry_num + 1
        if current_retry_num > max_retry_num:
            print(f'======= maximum retry num reached: {current_retry_num}')
            break

# event file
event_file = sys.argv[1]
with open(event_file, 'r') as f:
    events = [l.strip().split('|') for l in f.readlines() if not l.startswith('#')]

# data time range
time_before_origin = 500
time_after_origin = 1800

stationxml_dir = 'stations'
station_starttime = UTCDateTime('1970-01-01')
station_endtime = UTCDateTime('2030-01-01')

mseed_storage = 'waveforms'
log_dir = 'log'

domain = CircularDomain(latitude=54, longitude=4, minradius=0.0, maxradius=40.0)
# domain = CircularDomain(latitude=54, longitude=4, minradius=0.0, maxradius=10.0)
# domain = RectangularDomain(minlatitude=minlatitude, minlongitude=minlongitude, maxlatitude=maxlatitude, maxlongitude=maxlongitude)

for evt in events:
    print(f'\n{"#"*20} \n# Event: {evt}\n{"#"*20}')

    evnm = evt[0]
    evot = UTCDateTime(evt[1])
    evla = float(evt[3])
    evlo = float(evt[4])
    evdp = float(evt[5])
    evmag = float(evt[6])

    # evla = float(evt[3])
    # evlo = float(evt[4])
    # evdp = float(evt[5])
    # evmag = float(evt[6])

    data_t0 = evot-time_before_origin
    data_t1 = evot+time_after_origin

    restrictions = Restrictions(
        starttime=data_t0,
        endtime=data_t1,
        station_starttime=station_starttime,
        station_endtime=station_endtime,
        # chunklength_in_sec=300,
        reject_channels_with_gaps=False,
        minimum_length=0.95,
        minimum_interstation_distance_in_m=0,
        sanitize=False,
        # network="Z3", station="D159", location="00", channel="HHZ",
        # network="BQ", station="DREG", location="", channel="HHZ",
        # network="BQ", location="", # channel="HHZ",
        channel_priorities=["HH[ZNE12]", "BH[ZNE12]"],
        # location_priorities=["", "00", "10"]
        )

    event_log_dir = os.path.join(log_dir, evnm)
    if not os.path.exists(event_log_dir):
        os.makedirs(event_log_dir)

    download_data(domain, restrictions, mseed_storage, stationxml_dir, event_log_dir, max_retry_num=1, providers=['BGR'])
    # download_data(domain, restrictions, mseed_storage, stationxml_dir, event_log_dir, max_retry_num=2, providers=['BGR'])
    # download_data(domain, restrictions, mseed_storage, stationxml_dir, event_log_dir, max_retry_num=3, providers=['EIDA'])
