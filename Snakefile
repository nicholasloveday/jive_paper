rule all:
  input:
    "results/figures/dailypop.pdf",
    "results/station_maps/e_inland_stations.pdf",
    "results/figures/temperature.pdf",
    "results/station_maps/d_aus_temperature.pdf",
    "results/figures/wind_max_exp.pdf",
    "results/station_maps/a_wind_aus_stations.pdf",
    "results/figures/sw_wa_wind_bias.pdf",
    "results/station_maps/c_sw_wa_stations.pdf",
    "results/figures/wind_exposed_peaks.pdf",
    "results/station_maps/b_wind_exposed_peaks_stations.pdf"

rule fetch_data_from_zenodo:
  conda: "workflow/envs/get_data.yaml"
  log: "logs/fetch_data.log"
  output:
    "data_for_jive_paper.zip"
  shell: 
      "zenodo_get -r 11015211"

rule unzip_data_from_zenodo:
  input:
    "data_for_jive_paper.zip"
  output:
    "data/dailypop/AutoFcst_DailyPoP_12_20211201-20220228_inland.nc",
    "data/dailypop/obs_DailyPoP_20211201-20220228_inland.nc",
    "data/dailypop/Official_DailyPoP_00_20211201-20220228_inland.nc",
    "data/exposed_peaks/FCF_1_0_AutoFcst_WindMag_18_20171201-20180228_exposed_peaks.nc",
    "data/exposed_peaks/FCF_1_3_AutoFcst_WindMag_18_20171201-20180228_exposed_peaks.nc",
    "data/exposed_peaks/obs_WindMag_20171201-20180228_exposed_peaks.nc",
    "data/exposed_peaks/Official_WindMag_00_20171201-20180228_exposed_peaks.nc",
    "data/sw_wa_wind/AutoFcst_WindMag_12_20230601-20230830_sw_wa.nc",
    "data/sw_wa_wind/obs_WindMag_20230601-20230830sw_wa.nc",
    "data/sw_wa_wind/Official_WindMag_00_20230601-20230830_sw_wa.nc",
    "data/temperature/AutoFcst_MaxT_12_20200501-20200930_NT.nc",
    "data/temperature/obs_MaxT_20160501-20160930_NT.nc",
    "data/temperature/obs_MaxT_20191201-20200229.nc",
    "data/temperature/obs_MaxT_20200501-20200930_NT.nc",
    "data/temperature/Official_MaxT_00_20160501-20160930_NT.nc",
    "data/temperature/Official_MaxT_00_20191201-20200229.nc",
    "data/temperature/Official_MaxT_00_20200501-20200930_NT.nc",
    "data/temperature/PtOCF_MaxT_12_20160501-20160930_NT.nc",
    "data/temperature/PtOCF_MaxT_12_20191201-20200229.nc",
    "data/windmax_exp/AutoFcstMax_WindMag_12_20190901-20191130.nc",
    "data/windmax_exp/FCF_2_0_AutoFcst_WindMag_18_20190901-20191130.nc",
    "data/windmax_exp/obs_WindMagMaxInHour_20190901-20191130.nc",
    "data/windmax_exp/Official_WindMag_00_20190901-20191130.nc"
  shell:
    "unzip data_for_jive_paper.zip"


rule dailypop_calibration_figure:
  conda: "workflow/envs/notebooks.yaml"
  log: "logs/dailypop_calibration_figure.log"
  input: 
    "data/dailypop/Official_DailyPoP_00_20211201-20220228_inland.nc",
    "data/dailypop/AutoFcst_DailyPoP_12_20211201-20220228_inland.nc",
    "data/dailypop/obs_DailyPoP_20211201-20220228_inland.nc"
  output: 
    "results/figures/dailypop.pdf",
    "results/station_maps/e_inland_stations.pdf"
  notebook: "src/dailypop_calibration_fig.ipynb"

rule temperature_figure:
  conda: "workflow/envs/notebooks.yaml"
  log: "logs/temperature_figure.log"
  input:
    "data/climate/max_t_0.97.json",
    "data/temperature/AutoFcst_MaxT_12_20200501-20200930_NT.nc",
    "data/temperature/obs_MaxT_20160501-20160930_NT.nc",
    "data/temperature/obs_MaxT_20191201-20200229.nc",
    "data/temperature/obs_MaxT_20200501-20200930_NT.nc",
    "data/temperature/Official_MaxT_00_20160501-20160930_NT.nc",
    "data/temperature/Official_MaxT_00_20191201-20200229.nc",
    "data/temperature/Official_MaxT_00_20200501-20200930_NT.nc",
    "data/temperature/PtOCF_MaxT_12_20160501-20160930_NT.nc",
    "data/temperature/PtOCF_MaxT_12_20191201-20200229.nc"
  output: 
    "results/figures/temperature.pdf",
    "results/station_maps/d_aus_temperature.pdf"
  notebook: "src/temperature_fig.ipynb"

rule wind_max_exp_fig:
  conda: "workflow/envs/notebooks.yaml"
  log: "logs/wind_max_exp_fig.log"
  input: 
    "data/windmax_exp/AutoFcstMax_WindMag_12_20190901-20191130.nc",
    "data/windmax_exp/FCF_2_0_AutoFcst_WindMag_18_20190901-20191130.nc",
    "data/windmax_exp/obs_WindMagMaxInHour_20190901-20191130.nc",
    "data/windmax_exp/Official_WindMag_00_20190901-20191130.nc"
  output: 
    "results/figures/wind_max_exp.pdf",
    "results/station_maps/a_wind_aus_stations.pdf"
  notebook: "src/wind_max_exp_fig.ipynb"

rule wind_bias_figures:
  conda: "workflow/envs/notebooks.yaml"
  log: "logs/wind_bias_figures.log"
  input:
    "data/sw_wa_wind/AutoFcst_WindMag_12_20230601-20230830_sw_wa.nc",
    "data/sw_wa_wind/obs_WindMag_20230601-20230830sw_wa.nc",
    "data/sw_wa_wind/Official_WindMag_00_20230601-20230830_sw_wa.nc"
  output: 
    "results/figures/sw_wa_wind_bias.pdf",
    "results/station_maps/c_sw_wa_stations.pdf"
  notebook: "src/sw_wa_wind_fig.ipynb"


rule exposed_peaks_fig:
  conda: "workflow/envs/notebooks.yaml"
  log: "logs/exposed_peaks_fig.log"
  input:
    "data/exposed_peaks/FCF_1_0_AutoFcst_WindMag_18_20171201-20180228_exposed_peaks.nc",
    "data/exposed_peaks/FCF_1_3_AutoFcst_WindMag_18_20171201-20180228_exposed_peaks.nc",
    "data/exposed_peaks/obs_WindMag_20171201-20180228_exposed_peaks.nc",
    "data/exposed_peaks/Official_WindMag_00_20171201-20180228_exposed_peaks.nc"
  output: 
    "results/figures/wind_exposed_peaks.pdf",
    "results/station_maps/b_wind_exposed_peaks_stations.pdf"
  notebook: "src/exposed_peaks_fig.ipynb"
