#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Segy2netcdf: convert SEG-Y files to NetCDF files.
"""
import click
import obspy.io.segy.segy
#from obspy.io.segy.segy import _read_segy
import numpy as np
from netCDF4 import Dataset

@click.command()
@click.argument('segy_path', type=click.Path(exists=True, dir_okay=False))
@click.argument('netcdf_path', type=click.Path())
@click.option('--samples_dim_name', '-sdn', type=str,
        help='Name of trace samples dimension (usually Time or Depth)')
@click.option('-d', type=(str, int), multiple=True,
        help='Name and length (separated by a space) '
        'of other dimensions, in slowest to fastest order. If the name '
        'matches the name of a trace header field (using '
        'segyio.TraceField names), then the header values will be used '
        'for the dimension coordinates, otherwise the dimension '
        'coordinates will start at 0 and increment by 1. '
        'E.g. -d FieldRecord 10 -d GroupNumber 5, indicates that the '
        'data consists of 10 shot gathers, each with 5 receiver traces. '
        'As FieldRecord is a trace header name, the trace header values '
        'will be used as coordinates for that dimension. GroupNumber is '
        'not a header name, so the coordinates of that dimension will '
        'be 0 to 4.')
@click.option('--compress/--no-compress', default=False,
        help='turn on or off NetCDF compression (default off).')
@click.option('--verbose/--quiet', default=False,
        help='turn on or off verbose output (default off).')
#@click.option('--ignore-geometry', is_flag=True, default=False,
#        help='ignore geometry from the segy file (default False).')
def segy2netcdf_cli(segy_path, netcdf_path, samples_dim_name, d, compress, verbose):
  """Convert a SEG-Y file to a NetCDF file ()."""
  segy2netcdf(segy_path, netcdf_path, samples_dim_name, d, compress, verbose)

# global varibales
nsamples = 0
ntraces = 0
sample_interval = 0

segybin = obspy.io.segy.segy.SEGYBinaryFileHeader()
segy_binary_header_fields = [ (attr, type(getattr(segybin, attr)) ) for attr in dir(segybin) 
    if not callable(getattr(segybin, attr)) and not attr.startswith('__')]

#trace = obspy.io.segy.segy.SEGYTraceHeader()
#trace_header_fields = []
#trace_header_fields = [ (attr, type(getattr(trace, attr)) ) for attr in dir(trace) 
#    if not callable(getattr(trace, attr)) 
#    and not attr.startswith('__') # ]
#    and not attr.startswith('unassigned') ] # NOTE to avoid type conversion error, but this field might be usefull in some customized format.
trace_header_map = {
  'trace_sequence_number_within_line':                                              'TRACE_SEQUENCE_LINE'                   ,
  'trace_sequence_number_within_segy_file':                                         'TRACE_SEQUENCE_FILE'                   ,
  'original_field_record_number':                                                   'FieldRecord'                           ,
  'trace_number_within_the_original_field_record':                                  'TraceNumber'                           ,
  'energy_source_point_number':                                                     'EnergySourcePoint'                     ,
  'ensemble_number':                                                                'CDP'                                   ,
  'trace_number_within_the_ensemble':                                               'CDP_TRACE'                             ,
  'trace_identification_code':                                                      'TraceIdentificationCode'               ,
  'number_of_vertically_summed_traces_yielding_this_trace':                         'NSummedTraces'                         ,
  'number_of_horizontally_stacked_traces_yielding_this_trace':                      'NStackedTraces'                        ,
  'data_use':                                                                       'DataUse'                               ,
  'distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group':   'offset'                                ,
  'receiver_group_elevation':                                                       'ReceiverGroupElevation'                ,
  'surface_elevation_at_source':                                                    'SourceSurfaceElevation'                ,
  'source_depth_below_surface':                                                     'SourceDepth'                           ,
  'datum_elevation_at_receiver_group':                                              'ReceiverDatumElevation'                ,
  'datum_elevation_at_source':                                                      'SourceDatumElevation'                  ,
  'water_depth_at_source':                                                          'SourceWaterDepth'                      ,
  'water_depth_at_group':                                                           'GroupWaterDepth'                       ,
  'scalar_to_be_applied_to_all_elevations_and_depths':                              'ElevationScalar'                       ,
  'scalar_to_be_applied_to_all_coordinates':                                        'SourceGroupScalar'                     ,
  'source_coordinate_x':                                                            'SourceX'                               ,
  'source_coordinate_y':                                                            'SourceY'                               ,
  'group_coordinate_x':                                                             'GroupX'                                ,
  'group_coordinate_y':                                                             'GroupY'                                ,
  'coordinate_units':                                                               'CoordinateUnits'                       ,
  'weathering_velocity':                                                            'WeatheringVelocity'                    ,
  'subweathering_velocity':                                                         'SubWeatheringVelocity'                 ,
  'uphole_time_at_source_in_ms':                                                    'SourceUpholeTime'                      ,
  'uphole_time_at_group_in_ms':                                                     'GroupUpholeTime'                       ,
  'source_static_correction_in_ms':                                                 'SourceStaticCorrection'                ,
  'group_static_correction_in_ms':                                                  'GroupStaticCorrection'                 ,
  'total_static_applied_in_ms':                                                     'TotalStaticApplied'                    ,
  'lag_time_A':                                                                     'LagTimeA'                              ,
  'lag_time_B':                                                                     'LagTimeB'                              ,
  'delay_recording_time':                                                           'DelayRecordingTime'                    ,
  'mute_time_start_time_in_ms':                                                     'MuteTimeStart'                         ,
  'mute_time_end_time_in_ms':                                                       'MuteTimeEND'                           ,
  'number_of_samples_in_this_trace':                                                'TRACE_SAMPLE_COUNT'                    ,
  'sample_interval_in_ms_for_this_trace':                                           'TRACE_SAMPLE_INTERVAL'                 ,
  'gain_type_of_field_instruments':                                                 'GainType'                              ,
  'instrument_gain_constant':                                                       'InstrumentGainConstant'                ,
  'instrument_early_or_initial_gain':                                               'InstrumentInitialGain'                 ,
  'correlated':                                                                     'Correlated'                            ,
  'sweep_frequency_at_start':                                                       'SweepFrequencyStart'                   ,
  'sweep_frequency_at_end':                                                         'SweepFrequencyEnd'                     ,
  'sweep_length_in_ms':                                                             'SweepLength'                           ,
  'sweep_type':                                                                     'SweepType'                             ,
  'sweep_trace_taper_length_at_start_in_ms':                                        'SweepTraceTaperLengthStart'            ,
  'sweep_trace_taper_length_at_end_in_ms':                                          'SweepTraceTaperLengthEnd'              ,
  'taper_type':                                                                     'TaperType'                             ,
  'alias_filter_frequency':                                                         'AliasFilterFrequency'                  ,
  'alias_filter_slope':                                                             'AliasFilterSlope'                      ,
  'notch_filter_frequency':                                                         'NotchFilterFrequency'                  ,
  'notch_filter_slope':                                                             'NotchFilterSlope'                      ,
  'low_cut_frequency':                                                              'LowCutFrequency'                       ,
  'high_cut_frequency':                                                             'HighCutFrequency'                      ,
  'low_cut_slope':                                                                  'LowCutSlope'                           ,
  'high_cut_slope':                                                                 'HighCutSlope'                          ,
  'year_data_recorded':                                                             'YearDataRecorded'                      ,
  'day_of_year':                                                                    'DayOfYear'                             ,
  'hour_of_day':                                                                    'HourOfDay'                             ,
  'minute_of_hour':                                                                 'MinuteOfHour'                          ,
  'second_of_minute':                                                               'SecondOfMinute'                        ,
  'time_basis_code':                                                                'TimeBaseCode'                          ,
  'trace_weighting_factor':                                                         'TraceWeightingFactor'                  ,
  'geophone_group_number_of_roll_switch_position_one':                              'GeophoneGroupNumberRoll1'              ,
  'geophone_group_number_of_trace_number_one':                                      'GeophoneGroupNumberFirstTraceOrigField',
  'geophone_group_number_of_last_trace':                                            'GeophoneGroupNumberLastTraceOrigField' ,
  'gap_size':                                                                       'GapSize'                               ,
  'over_travel_associated_with_taper':                                              'OverTravel'                            ,
  'x_coordinate_of_ensemble_position_of_this_trace':                                'CDP_X'                                 ,
  'y_coordinate_of_ensemble_position_of_this_trace':                                'CDP_Y'                                 ,
  'for_3d_poststack_data_this_field_is_for_in_line_number':                         'INLINE_3D'                             ,
  'for_3d_poststack_data_this_field_is_for_cross_line_number':                      'CROSSLINE_3D'                          ,
  'shotpoint_number':                                                               'ShotPoint'                             ,
  'scalar_to_be_applied_to_the_shotpoint_number':                                   'ShotPointScalar'                       ,
  'trace_value_measurement_unit':                                                   'TraceValueMeasurementUnit'             ,
  'transduction_constant_mantissa':                                                 'TransductionConstantMantissa'          ,
  'transduction_constant_exponent':                                                 'TransductionConstantPower'             ,
  'transduction_units':                                                             'TransductionUnit'                      ,
  'device_trace_identifier':                                                        'TraceIdentifier'                       ,
  'scalar_to_be_applied_to_times':                                                  'ScalarTraceHeader'                     ,
  'source_type_orientation':                                                        'SourceType'                            ,
  'source_energy_direction_mantissa':                                               'SourceEnergyDirectionMantissa'         ,
  'source_energy_direction_exponent':                                               'SourceEnergyDirectionExponent'         ,
  'source_measurement_mantissa':                                                    'SourceMeasurementMantissa'             ,
  'source_measurement_exponent':                                                    'SourceMeasurementExponent'             ,
  'source_measurement_unit':                                                        'SourceMeasurementUnit'                 ,
#  'unassigned':                                                                     'UnassignedInt'                         ,
}

def segy2netcdf(segy_path, netcdf_path, samples_dim_name=None, d=(),
        compress=False, verbose=False,
        #        ignore_geometry=False,
        ):
  """Convert a SEG-Y file to a NetCDF file.

  Args:
    segy_path: A string specifying the path to input SEG-Y file.
    netcdf_path: A string specifying the path to output NetCDF file
    samples_dim_name: An optional string specifying the name of trace
      samples dimension (usually Time or Depth). Default SampleNumber.
    d: An optional tuple of tuples of the form (string, int), where the
      string specifies a dimension name, and the int specifies the number
      of entries in that dimension. There should be one of these inner
      tuples for each dimension, excluding the trace sample dimension,
      in slowest to fastest order. Any dimensions not accounted for,
      either because d was not specified, or it was specified but did
      not account for all of the traces in the SEG-Y file, a Traces
      dimension will be used for the remainder.
    compress: An optional boolean flag indicating whether NetCDF
      compression should be used. Default False.
    verbose: An optional boolean flag indicating whether to print
      progress. Default False.
  """

  global nsamples
  global sample_interval
  global ntraces

  # set default name for trace samples dimension
  if not samples_dim_name:
    samples_dim_name = 'SampleNumber'

  segy = obspy.io.segy.segy._read_segy(segy_path)

  nsamples = segy.binary_file_header.number_of_samples_per_data_trace
  sample_interval = segy.binary_file_header.sample_interval_in_microseconds
  ntraces = len(segy.traces)

  dim_names, dim_lens = _make_dim_name_len(segy, samples_dim_name, nsamples, d)
  dims_ntraces = _count_traces_in_user_dims(d)

  _check_user_dims(dims_ntraces, ntraces)

  _fill_missing_dims(dims_ntraces, ntraces, dim_names, dim_lens)

  rootgrp = Dataset(netcdf_path, 'w', format='NETCDF4')
  _create_dimensions(dim_names, dim_lens, rootgrp)
  variables = _create_variables(rootgrp, dim_names, compress)
  _set_attributes(segy, rootgrp)
  _copy_data(segy, variables, dim_names, dim_lens, verbose)

  rootgrp.close()


def _make_dim_name_len(segy, samples_dim_name, nsamples, d):
  """Make dim_names and dim_lens lists of dimension names and lens.

  Args:
    samples_dim_name: A string specifying the name to use for the trace
      samples dimension. Usually 'Time' or 'Depth'
    ns: An int specifying the number of samples per trace
    d: A tuple of tuples of the form (string, int), where the string
      specifies a dimension name, and the int specifies the number of
      entries in that dimension.

  Returns:
    dim_names: A list with the dimension names
    dim_lens: A list with the lengths of the dimensions
  """
  global trace_header_map
  header_names = list(trace_header_map.values())
  header_fields = list(trace_header_map.keys())

  dim_names = []
  dim_lens = []
  dim_unique_lens = []
  for dim in d:
    unique_len = -1
    if dim[0] in header_names:
      field = header_fields[header_names.index(dim[0])]
      unique_len = np.unique(np.array([ getattr(tr.header, field) for tr in segy.traces ])).size
    dim_unique_lens.append(unique_len)
    dim_names.append(dim[0])
    dim_lens.append(dim[1])
  dim_names.append(samples_dim_name)
  dim_lens.append(nsamples)
  print("dim_names: ", dim_names)
  print("dim_lens: ", dim_lens)
  print("dim_unique_lens(without_sample_dim): ", dim_unique_lens)
  return dim_names, dim_lens


def _count_traces_in_user_dims(d):
  """Determine how many traces are accounted for by the dimensions provided
     by the user.
  """
  dims_ntraces = 1
  for dim in d:
    dims_ntraces *= dim[1]
  return dims_ntraces


def _check_user_dims(dims_ntraces, ntraces):
  """Check that the lengths of the user-provided dimensions make sense."""
  if dims_ntraces > ntraces:
    raise ValueError('supplied dimensions imply {} traces, '
             'but only {} in file'.format(dims_ntraces, ntraces))
  if dims_ntraces < 0:
    raise ValueError('supplied dimensions imply {} traces'
             .format(dims_ntraces))
  if (dims_ntraces == 0) and (ntraces > 0):
    raise ValueError('supplied dimensions imply {} traces, '
             'but {} in file'.format(dims_ntraces, ntraces))
  if (dims_ntraces != 0) and (ntraces % dims_ntraces != 0):
    raise ValueError('supplied dimensions imply {} traces, '
             'but this does not divide into the {} traces in the '
             'file'.format(dims_ntraces, ntraces))


def _fill_missing_dims(dims_ntraces, ntraces, dim_names, dim_lens):
  """Make a new dimension (if necessary) for unaccounted traces."""
  if dims_ntraces < ntraces:
    dim_names.insert(0, 'Traces')
    dim_lens.insert(0, int(ntraces / dims_ntraces))


def _create_dimensions(dim_names, dim_lens, rootgrp):
  """Create the dimensions in the NetCDF file."""
  for dim in zip(dim_names, dim_lens):
    rootgrp.createDimension(dim[0], dim[1])


def _create_variables(rootgrp, dim_names, compress):
  """Create variables in the NetCDF file.

     The trace data, Time/Depth dimension, and trace headers, are all
     created as variables.
  """
  variables = []
  # Trace data
  variables.append(rootgrp.createVariable('Samples', 'f4', tuple(dim_names),
                      zlib=compress))
  # Time/Depth dimension
  variables.append(rootgrp.createVariable(dim_names[-1], 'f4',
                      dim_names[-1], zlib=compress))
  # Other dimensions
  variables += _create_traceheader_variables(rootgrp, dim_names, compress)
  return variables


def _create_traceheader_variables(rootgrp, dim_names, compress):
  """Create NetCDF variables for each trace header field.

     Fields that are used as dimensions are only the length of that
     dimension, others have one entry for every trace.
  """
  #global trace_header_fields
  global trace_header_map

  fieldtype = 'i4'
  variables = []
  #for fieldname, fieldtype in trace_header_fields:
  for fieldname in trace_header_map:
    # map fieldname to varname if exists
    #if fieldname in trace_header_map: 
    #  varname = trace_header_map[fieldname]
    #else:
    #  varname = fieldname
    varname = trace_header_map[fieldname]
    # for variables that are dimensions of the dataset, they should be the
    # size of their dimension. All others should be the size of the dataset
    # (excluding the trace samples dimension)
    if varname in dim_names:
      v = rootgrp.createVariable(varname, fieldtype, varname, zlib=compress)
      v.fieldname = fieldname
      variables.append(v)
    else:
      v = rootgrp.createVariable(varname, fieldtype, tuple(dim_names[:-1]), zlib=compress)
      v.fieldname = fieldname
      variables.append(v)

  return variables


def _set_attributes(segy, rootgrp):
  """Copy the file headers (binary and text) to the NetCDF file."""
  global segy_binary_header_fields

  for fieldname, fieldtype in segy_binary_header_fields:
    field_value = getattr(segy.binary_file_header, fieldname)
    setattr(rootgrp, fieldname, field_value)
  rootgrp.textual_file_header = segy.textual_file_header
  #if segy.ext_headers:
  #  rootgrp.ext_headers = segy.text[1]

def _copy_data(segy, variables, dim_names, dim_lens, verbose):
  """Copy the data to the NetCDF file.

     Trace data, Time/Depth dimension indices, and trace header values are
     copied.

  """
  global nsamples
  global sample_interval

  n_trace_dims = len(dim_names[:-1])

  for v in variables:
    if v.name == 'Samples':
      if verbose: click.echo('copying trace data')
      v[:] = np.array([ tr.data for tr in segy.traces ])
    elif v.name == dim_names[-1]:
      if verbose: click.echo('copying time/depth indices')
      v[:] = np.arange(0, nsamples) * sample_interval
    else:
      if verbose: click.echo('copying {}'.format(v.name))
      data = np.reshape( np.array([ getattr(tr.header, v.fieldname) for tr in segy.traces]), dim_lens[:-1])
      # Headers used as dimensions will only copy from traces that should contain unique values for them, while other headers will copy from every trace.
      dims = [0,] * n_trace_dims
      for d in v.dimensions:
        d_idx = dim_names.index(d)
        dims[d_idx] = slice(None)
      try:
        v[:] = data[tuple(dims)]
      except:
        print("varaible: ", v, v.datatype) 
        print("input data: ", data.shape, data.dtype )
        raise

if __name__ == '__main__': 
  segy2netcdf_cli()
