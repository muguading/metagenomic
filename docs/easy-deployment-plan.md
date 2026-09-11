# Easy Deployment Plan

本方案的目标不是把整个病原分析体系伪装成一个“下一步下一步”的小软件。那是漂亮的谎言。这个项目真正难部署的地方不在 Flask，而在 Conda 环境、第三方生信工具、参考数据库和大体积资产。所以部署必须分层，把轻的做成一键，把重的做成可检查、可挂载、可替换。

## 结论

默认交付形态采用三层结构：

1. **Portal 应用层**：Flask + gunicorn + SQLite，Docker Compose 或 systemd 一键启动。
2. **运行环境层**：Conda/Mamba 环境独立安装，由管理页扫描并绑定环境名。
3. **数据资产层**：数据库、参考序列、第三方工具包不进镜像，统一挂载到固定目录并由部署清单检查。

这条路最朴素，也最不容易翻车。把数据库塞进镜像会得到一个巨大、难更新、难分发、难排错的怪物；把所有工具都要求用户手工配置，又会把软件交付退化成口口相传。

## 三个方案

### 方案 A：渐进方案，先把 Portal 部署简单化

适合：院内服务器、单机工作站、演示环境、第一版对外交付。

交付内容：

- `deployment/Dockerfile.portal`
- `deployment/docker-compose.portal.yml`
- `deployment/portal.env.example`
- `deploy_bac_analysis_portal_ubuntu.sh` 作为 systemd/nginx 备用路径继续保留

部署体验：

```bash
cp deployment/portal.env.example deployment/portal.env
python3.10 -c 'import secrets; print(secrets.token_urlsafe(48))'
python3.10 scripts/check_deployment.py --env-file deployment/portal.env
docker compose -f deployment/docker-compose.portal.yml --env-file deployment/portal.env up -d --build
```

用户随后访问 `http://服务器IP:5055`，进入管理页选择：

- 部署基准目录
- 主分析脚本
- Conda 安装路径
- 数据库部署目录
- 各流程 Conda 环境名

优点：

- 改动小，最快落地。
- Portal 和数据库资产解耦，升级应用不需要重打几十 GB 数据。
- 失败边界清楚：Web 起不来、Conda 找不到、数据库缺失分别排查。

缺点：

- 首次安装 Conda 环境和数据库仍需要管理员参与。
- Docker 容器调用宿主 Conda 工具链时，需要谨慎规划挂载路径和权限；生产场景更推荐把 Portal 跑在服务器本机 systemd 下。

判断：这是现在应该做的默认方案。别为了所谓“全自动”牺牲可维护性。

### 方案 B：激进方案，做离线安装器

适合：无公网医院、客户现场、需要 U 盘交付的 Windows/Linux 工作站。

交付内容：

- `pathogen-workbench-offline-YYYYMMDD.tar.zst`
- 内含源码、wheelhouse、conda env yml、最小数据库、校验清单、安装脚本
- `install.sh` / `install.ps1` 负责安装 Python 环境、恢复配置、注册服务

部署体验：

```bash
tar --zstd -xf pathogen-workbench-offline-YYYYMMDD.tar.zst
cd pathogen-workbench-offline
bash install.sh --prefix /opt/pathogen-workbench --database-root /data/pathogen-db
```

优点：

- 对客户最像“软件安装”。
- 不依赖公网，交付确定性强。
- 可以做版本冻结和验收清单。

缺点：

- 制作成本明显更高。
- 每次数据库或工具链更新都要重新生成资产包。
- 安装器必须做严格校验，否则只是把手工错误换成自动错误。

判断：第二阶段做。等 Portal 默认部署稳定后，再把它封装成离线安装器。

### 方案 C：理想方案，平台化部署

适合：多科室、多用户、多队列、长期运维。

交付内容：

- Portal 服务容器化
- 任务执行器独立成 worker
- 数据库资产进入只读共享卷
- PostgreSQL 替代 SQLite
- MinIO/NAS 管理输入输出
- Prometheus/Grafana 做资源与任务监控

优点：

- 多用户、多任务、可观测性和扩展性最好。
- 应用升级、任务执行和数据资产各自独立。

缺点：

- 这是平台，不是“小软件”。
- 运维门槛高，必须有人负责服务器、存储、权限、备份和监控。

判断：预算、团队和使用规模都上来之后再做。现在硬上，只会把交付周期拖进泥里。

## 推荐目录

```text
/opt/pathogen-workbench/          # 应用代码
/opt/pathogen-workbench/.venv_web # Python Web 环境，systemd 方案使用
/data/pathogen-workbench/state    # SQLite、缓存、运行状态
/data/pathogen-workbench/tasks    # 任务输入、日志、输出
/data/pathogen-db                 # 数据库与参考资产
/opt/miniconda3                   # Conda/Mamba 安装
```

## 执行权限与数据目录

生产安装器默认创建不可登录的 `pathogen-workbench` 系统账号。Portal 及其分析子进程均以该账号运行；不要使用 `root`、个人登录账号或共享管理员账号启动服务。

- `/opt/pathogen-workbench/app`、Conda 运行时和 `/data/pathogen-db` 应由管理员维护，并对服务账号只读。
- `/data/pathogen-workbench/state` 与 `/data/pathogen-workbench/tasks` 是服务账号唯一默认可写目录，权限为 `0750`。
- 分析输出必须置于由该账号管理的批准输出根目录；不要把个人主目录、系统目录或参考库目录作为输出位置。
- systemd 使用 `NoNewPrivileges`、私有 `/tmp`、只读系统/代码/参考库路径和最小 `ReadWritePaths`。如新增数据盘，必须显式加入服务的可写路径，并按同一账号/权限规则创建。

Docker Compose 同样使用 UID/GID `10001` 的非 root 账号、只读根文件系统、临时 `/tmp`、只读参考库挂载及 `no-new-privileges`。宿主机的 state/task 挂载目录应事先 `chown 10001:10001`。

Docker Compose 方案中，源码目录挂载为 `/app`，状态目录挂载为 `/data/pathogen-workbench/state`，数据库目录挂载为 `/data/pathogen-db`。

## 推荐部署命令

### 完整离线包（推荐给真实分析）

适合目标服务器不能联网、但需要部署后直接跑真实样本的生产交付。先在 Ubuntu 构建服务器上准备 Conda 环境、数据库和工具资产，再生成一个可拷贝到移动硬盘的单目录发布包。

构建服务器：

```bash
bash scripts/build_linux_full_bundle.sh \
  --version 2026.06 \
  --source-root /path/to/metagenomic \
  --conda-root /opt/miniconda3 \
  --database-root /data/pathogen-db \
  --output /mnt/usb/pathogen-workbench-full-linux \
  --archive
```

发布目录结构：

```text
pathogen-workbench-full-linux/
├── install.sh
├── manifest.json
├── app/
├── conda-runtime/
├── conda-packs/
├── database/
├── soft/
├── public/
└── smoke-tests/
```

目标 Ubuntu 服务器：

```bash
cd /mnt/usb
sudo bash pathogen-workbench-full-linux/install.sh --yes
```

目标服务器基础系统需要预装 `python3`、`tar`、`sha256sum`、`zstd`、`systemd` 和 `nginx`。这是离线安装的底线：安装器不在目标服务器联网 `apt install`。

安装器会校验 manifest、复制应用和数据库、解压 `conda-pack` 环境、运行 `conda-unpack`、写入 `portal.env`、初始化 Portal 部署设置、注册 systemd/nginx，并执行：

```bash
/opt/pathogen-workbench/conda/envs/web_runtime/bin/python \
  /opt/pathogen-workbench/app/scripts/check_deployment.py \
  --project-root /opt/pathogen-workbench/app \
  --env-file /opt/pathogen-workbench/app/deployment/portal.env \
  --profile full-analysis
```

这条路线不在目标服务器执行 `conda env create`，也不要求目标服务器访问公网。

### Docker Compose

适合先跑通 Web 工作台、演示、轻量服务器部署。

```bash
cp deployment/portal.env.example deployment/portal.env
python3.10 -c 'import secrets; print(secrets.token_urlsafe(48))'
```

把生成的随机值写入 `PORTAL_SECRET_KEY`，把 `PORTAL_INITIAL_ADMIN_PASSWORD` 改成强密码，然后执行：

```bash
python3.10 scripts/check_deployment.py --env-file deployment/portal.env
docker compose -f deployment/docker-compose.portal.yml --env-file deployment/portal.env up -d --build
```

### Ubuntu systemd + nginx

适合生产服务器，尤其是真实分析需要直接访问宿主机 Conda 和数据库资产时。

```bash
cp deployment/portal.env.example deployment/portal.env
python3.10 -c 'import secrets; print(secrets.token_urlsafe(48))'
bash deploy_bac_analysis_portal_ubuntu.sh
```

脚本会：

- 创建 `.venv_web`
- 安装 Web 依赖和 gunicorn
- 生成 `run_analysis_portal_linux.sh`
- 生成 systemd service 模板
- 生成 nginx site 模板
- 让 service 读取同一份 `deployment/portal.env`

安装 service 前先跑：

```bash
./.venv_web/bin/python scripts/check_deployment.py --env-file deployment/portal.env
```

## 环境变量

生产环境必须显式配置：

```bash
PORTAL_MODE=production
PORTAL_SECRET_KEY=<随机强密钥>
PORTAL_INITIAL_ADMIN_PASSWORD=<首次管理员强密码>
PORTAL_COOKIE_SECURE=0
PORTAL_DB_PATH=/data/pathogen-workbench/state/bac_analysis_portal.sqlite3
BAC_ANALYSIS_TASK_ROOT=/data/pathogen-workbench/tasks
APP_HOST=127.0.0.1
APP_PORT=5055
GUNICORN_WORKERS=2
```

如果前面有 HTTPS 反向代理，`PORTAL_COOKIE_SECURE` 应设为 `1`。如果只是内网 IP + 端口测试，先设为 `0`，别假装自己有 HTTPS。

真实分析常用路径：

```bash
META_DATABASE_ROOT=/data/pathogen-db
META_KRAKEN_DB=/data/pathogen-db/kraken2
META_VIRUS_KRAKEN_DB=/data/pathogen-db/virus_kraken2
META_VIRSORTER2_DB=/data/pathogen-db/virsorter2
META_CHECKV_DB=/data/pathogen-db/checkv-db
META_GENOMAD_DB=/data/pathogen-db/genomad
META_MOBILEOG_DB=/data/pathogen-db/mobileOG-db
META_MOBILEOG_META=/data/pathogen-db/mobileOG-db-beatrix.csv
```

## 最小验收清单

Portal 验收：

- 登录页能打开。
- 生产模式下不能使用 `admin / admin123`。
- 管理页能保存部署基准目录、数据库目录和 Conda 环境。
- 创建一个演示任务后，任务状态、日志和报告页能正常展示。

分析环境验收：

- `conda` 能被系统检测到。
- 所有必填环境名存在。
- 数据库目录可读。
- 每个工作站至少有一个小样本 smoke test。
- 失败日志能在 Portal 中看到，而不是只在服务器终端里消失。

可先用体检脚本做机器级检查：

```bash
python3.10 scripts/check_deployment.py --env-file deployment/portal.env --strict-assets
python3.10 scripts/check_deployment.py --env-file deployment/portal.env --check-default-conda-envs
python3.10 scripts/check_deployment.py --env-file deployment/portal.env --json
```

升级验收：

- 升级前备份 SQLite 和 `analysis_tasks/`。
- 新版本启动后旧任务列表仍可访问。
- 管理页配置不丢失。
- 离线更新源必须包含有效部署目录结构。

## 下一步

1. 用 `scripts/check_deployment.py` 固化每次交付前的验收输出。
2. 在 Portal 管理页增加“部署体检”按钮，把 `--profile full-analysis --json` 的检查结果可视化。
3. 把真实分析的最小 smoke test 数据集纳入交付包，避免只验证 Web 没验证流程。

真正值得追求的不是“安装过程看起来短”，而是失败时用户知道问题在哪里。部署体验的高级感，不是藏住复杂性，而是把复杂性驯服到一个清晰、可恢复的流程里。
