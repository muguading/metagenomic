# Batch Import Precheck Test Data

这组文件用于测试“批量导入预检 -> 正式导入”的闭环。

## 样本数据库

- `samples_valid.csv`: 应通过预检。正式导入会写入 2 条测试样本记录。
- `samples_invalid.csv`: 应无法通过预检，用于测试缺字段、空值和不存在路径的错误提示。
- `samples_with_a_very_long_filename_for_card_overflow_regression_testing_2026_05_19_valid.csv`: 内容与成功样本一致，用于回归测试长文件名不会撑破卡片。

## 参考基因组数据库

- `reference_host_valid.csv`: 应通过宿主参考基因组预检。正式导入会写入 1 条宿主测试记录。
- `reference_pathogen_valid.csv`: 应通过病原参考基因组预检。正式导入会写入 1 条病原测试记录。
- `reference_invalid.csv`: 应无法通过预检，用于测试缺少物种名、基因组名和不存在 FASTA 路径的错误提示。

## 注意

成功文件如果点击“开始批量导入”，会向当前数据库写入测试记录。只想测试预检和按钮状态时，不要点击正式导入。
