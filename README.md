# SGBM-SIMD

## 运行
./build/example/Release/sgbm.exe ./data/case1/imL.png ./data/case1/imR.png

## result
```
Loading Views...Done!
w = 1000, h = 750, d = [0,64]

SGM Initializing...
SGM Initializing Done! Timing : 0.236000 s

SGM Matching...
computing cost! timing :        0.076000 s
sgm_util::CostAggregateDagonal_1 Done! Timing : 0.271000 s

sgm_util::CostAggregateDagonal_2 Done! Timing : 0.284000 s

sgm_util::CostAggregateUpDown Done! Timing : 0.317000 s

sgm_util::CostAggregateDagonal_1 Done! Timing : 0.321000 s

sgm_util::CostAggregateDagonal_2 Done! Timing : 0.345000 s

sgm_util::CostAggregateLeftRight Done! Timing : 0.371000 s

sgm_util::CostAggregateUpDown Done! Timing : 0.372000 s

sgm_util::CostAggregateLeftRight Done! Timing : 0.402000 s

cost aggregating! timing :      0.479000 s
computing disparities! timing : 0.111000 s
postprocessing! timing :        0.206000 s

SGM Matching...Done! Timing :   0.876000 s
```