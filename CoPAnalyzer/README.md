# CoP Analyzer（CoP指標計算プログラム）
This program calculate stabilometry indices from CoP time series.

## Prerequisite
- OpenSim 4.0
- CMake (version >= 3.2)

## How to use
1. Put the CoP time series `jiang_CoPpositon_pos.sto` into the directory.
1. Run CMake with the `CMakeLists.txt`.
1. Build the project.
1. Execute the artifact `CoPAnalyzer`.
1. `results.csv` will be output.

## Settings
- `DT` (l.157) is the sampling period.
- `sto.crop(start,end)` (l.161) sets the time period.

## References
1. [今岡薫; 村瀬仁; 福原美穂. 重心動揺検査における健常者データの集計. Equilibrium research, 1997, 56.12Supplement: 1-84.](https://www.jstage.jst.go.jp/article/jser1971/56/12Supplement/56_12Supplement_1/_article/-char/ja)
