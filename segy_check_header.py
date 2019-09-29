#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Segy2netcdf: convert SEG-Y files to NetCDF files.
"""
import click
import pprint
import obspy.io.segy.segy
import numpy as np

@click.command()
@click.argument('segy_path', type=click.Path(exists=True, dir_okay=False))
@click.option('-d', type=str, multiple=True,
    default=('CDP_X', 'CDP_Y', 'FieldRecord', 'TraceNumber'),
    help='Name of  dimensions. If the name '
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
def segy_check_header_click(segy_path, d):
  """Convert a SEG-Y file to a NetCDF file ()."""
  segy_check_header(segy_path, d)

# global varibales
nsamples = 0
ntraces = 0
sample_interval = 0

#segybin = obspy.io.segy.segy.SEGYBinaryFileHeader()
#segy_binary_header_fields = [ (attr, type(getattr(segybin, attr)) ) for attr in dir(segybin) 
#    if not callable(getattr(segybin, attr)) and not attr.startswith('__')]

trace = obspy.io.segy.segy.SEGYTraceHeader()
trace_header_fields = [ (attr, type(getattr(trace, attr)) ) for attr in dir(trace) 
    if not callable(getattr(trace, attr)) 
    and not attr.startswith('__') # ]
    and not attr.startswith('unassigned') ] # NOTE to avoid type conversion error, but this field might be usefull in some customized format.

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
  'unassigned':                                                                     'UnassignedInt'                         ,
}

def segy_check_header(segy_path, d=()):
  """Convert a SEG-Y file to a NetCDF file.

  Args:
    segy_path: A string specifying the path to input SEG-Y file.
    d: An optional tuple of tuples of the form (string, int), where the
      string specifies a dimension name, and the int specifies the number
      of entries in that dimension. There should be one of these inner
      tuples for each dimension, excluding the trace sample dimension,
      in slowest to fastest order. Any dimensions not accounted for,
      either because d was not specified, or it was specified but did
      not account for all of the traces in the SEG-Y file, a Traces
      dimension will be used for the remainder.
  """
  global trace_header_map
  header_names = list(trace_header_map.values())
  header_fields = list(trace_header_map.keys())

  # set default name for trace samples dimension
  segy = obspy.io.segy.segy._read_segy(segy_path, headonly=True)
  #segy = obspy.io.segy.segy.iread_segy(segy_path, headonly=True)

  print(str(segy.binary_file_header))

  print("\n====== TEXTUAL_FILE_HEADER ======")
  pprint.pprint(segy.textual_file_header, width=84)

  print("\n====== DIMENSIONS ======")
  nsamples = segy.binary_file_header.number_of_samples_per_data_trace
  sample_interval = segy.binary_file_header.sample_interval_in_microseconds
  ntraces = len(segy.traces)
  print("number_of_samples_per_data_trace: ", nsamples)
  print("sample_interval_in_microseconds: ", sample_interval)
  print("number_of_data_trace: ", ntraces)

  dim_names = []
  dim_unique_lens = []
  for dim in d:
    unique_len = -1
    if dim in header_names:
      field = header_fields[header_names.index(dim)]
      data = np.array([ getattr(tr.header, field) for tr in segy.traces ])
      print(dim, '[0:5]= ', data[0:10])
      unique_len = np.unique(data).size
    dim_unique_lens.append(unique_len)
    dim_names.append(dim)
  print("dim_names: ", dim_names)
  print("dim_unique_lens: ", dim_unique_lens)


if __name__ == '__main__':
  segy_check_header_click()
