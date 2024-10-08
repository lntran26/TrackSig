# env: windows: micromamba activate tracksig_analysis

python simulation_analysis.py S8_S28_1 0.2 sum
python simulation_analysis.py S8_S28_1 0.2 bin 30
python simulation_analysis.py S2_S26_2 0.3 sum
python simulation_analysis.py S2_S26_2 0.4 bin 30
python simulation_analysis.py S2_S26_2 0.4 bin 49 # last bin due to 0-index
python simulation_analysis.py S2_S26_2 0.025 range 15 33

# ran after bug fix
python simulation_analysis.py S12_S22_3 0.2 sum
python simulation_analysis.py S12_S22_3 0.2 bin 27
