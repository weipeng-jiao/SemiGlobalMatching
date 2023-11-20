# SGBM Cmake 编译工程

## 运行
./build/example/Release/sgbm.exe ./data/case1/imL.png ./data/case1/imR.png

## result
```
Loading Views...Done!
w = 1000, h = 750, d = [0,64]

SGM Initializing...
SGM Initializing Done! Timing : 0.095000 s

SGM Matching...
computing cost! timing :        0.413000 s
cost aggregating! timing :      1.485000 s
computing disparities! timing : 0.085000 s
postprocessing! timing :        0.178000 s

SGM Matching...Done! Timing :   2.164000 s
```