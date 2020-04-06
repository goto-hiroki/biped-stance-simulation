# Objective closed stance optimization（閉脚立位の目標姿勢最適化）
This program optimize objective posture of closed stance on OpenSim. このプログラムは，OpenSim上の人体モデルにおいて，閉脚（閉足）立位の目標姿勢を最適化します．閉脚立位とは，両足のつま先同士とかかと同士をくっつけて立つ姿勢です．

## Prerequisite
- OpenSim 4.0
- CMake (version >= 3.2)
- [model data](https://figshare.com/articles/Source_code/7706903)[2]

## How to use
1. Put the model data into the directory.
2. Run CMake with the `CMakeLists.txt`.
3. Build the project.
4. Execute the artifact `ObjectivePostureOptimization`.
5. Objective posture `{minimum of objective function}.sto` will be output.

## References
1. [Jiang P, Chiba R, Takakusaki K, Ota J (2016) Generation of the Human Biped Stance by a Neural Controller Able to Compensate Neurological Time Delay. PLoS ONE 11(9): e0163212.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0163212)
2. [Kaminishi K, Jiang P, Chiba R, Takakusaki K, Ota J (2019) Postural control of a musculoskeletal model against multidirectional support surface translations. PLoS ONE 14(3): e0212613.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0212613#references)
