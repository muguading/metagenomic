#!/bin/bash
#
# cgMLST.org Panel 下载脚本 (Shell 版本)
# 使用 curl 进行下载
#

set -e

# 配置
OUTPUT_DIR="cgmlst_panels"
BASE_URL="https://www.cgmlst.org/ncs/api"
DELAY=2

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "cgMLST.org Panel 下载工具"
echo "========================================"
echo "输出目录: $OUTPUT_DIR"
echo ""

# 函数：安全的文件名
safe_filename() {
    echo "$1" | sed 's/[\\/<>|:"*? ]/_/g' | sed 's/__*/_/g' | sed 's/^_//;s/_$//'
}

# 函数：获取 schemes 列表
get_schemes() {
    echo "正在获取 schemes 列表..."
    curl -s -L \
        -H "User-Agent: Mozilla/5.0" \
        -H "Accept: application/json" \
        "$BASE_URL/schemes/"
}

# 函数：下载 scheme 详情
download_scheme_detail() {
    local scheme_id=$1
    local output_dir=$2
    curl -s -L \
        -H "User-Agent: Mozilla/5.0" \
        -H "Accept: application/json" \
        "$BASE_URL/scheme/$scheme_id" \
        -o "$output_dir/scheme_details.json"
}

# 函数：下载 alleles
download_alleles() {
    local scheme_id=$1
    local output_dir=$2

    local endpoints=(
        "$BASE_URL/scheme/$scheme_id/alleles/fasta"
        "$BASE_URL/scheme/$scheme_id/fasta"
    )

    for url in "${endpoints[@]}"; do
        echo "  尝试 alleles: $url"
        if curl -s -L -H "User-Agent: Mozilla/5.0" "$url" -o "$output_dir/alleles.fasta.tmp"; then
            # 检查文件是否有效（至少包含一个 > 字符）
            if head -1 "$output_dir/alleles.fasta.tmp" | grep -q ">"; then
                mv "$output_dir/alleles.fasta.tmp" "$output_dir/alleles.fasta"
                echo "    ${GREEN}✓ 下载成功${NC}"
                return 0
            fi
        fi
        rm -f "$output_dir/alleles.fasta.tmp"
        sleep 1
    done
    return 1
}

# 函数：下载 profiles
download_profiles() {
    local scheme_id=$1
    local output_dir=$2

    local endpoints=(
        "$BASE_URL/scheme/$scheme_id/profiles"
        "$BASE_URL/scheme/$scheme_id/profiles.tsv"
    )

    for url in "${endpoints[@]}"; do
        echo "  尝试 profiles: $url"
        if curl -s -L -H "User-Agent: Mozilla/5.0" "$url" -o "$output_dir/profiles.tsv.tmp"; then
            # 检查文件是否有内容
            if [ -s "$output_dir/profiles.tsv.tmp" ] && [ $(wc -l < "$output_dir/profiles.tsv.tmp") -gt 1 ]; then
                mv "$output_dir/profiles.tsv.tmp" "$output_dir/profiles.tsv"
                local lines=$(wc -l < "$output_dir/profiles.tsv")
                echo "    ${GREEN}✓ 下载成功 ($lines 行)${NC}"
                return 0
            fi
        fi
        rm -f "$output_dir/profiles.tsv.tmp"
        sleep 1
    done
    return 1
}

# 函数：处理单个 scheme
process_scheme() {
    local scheme_id=$1
    local scheme_name=$2
    local species=$3
    local index=$4
    local total=$5

    # 创建安全的目录名
    local safe_name=$(safe_filename "${scheme_id}_${species}_${scheme_name}")
    local scheme_dir="$OUTPUT_DIR/$safe_name"

    echo "[$index/$total] $species - $scheme_name (ID: $scheme_id)"

    # 检查是否已存在
    if [ -f "$scheme_dir/.completed" ]; then
        echo "  ${YELLOW}已存在，跳过${NC}"
        return 0
    fi

    mkdir -p "$scheme_dir"

    local success=false

    # 下载详情
    if download_scheme_detail "$scheme_id" "$scheme_dir"; then
        success=true
    fi

    # 下载 alleles
    if download_alleles "$scheme_id" "$scheme_dir"; then
        success=true
    fi

    # 下载 profiles
    if download_profiles "$scheme_id" "$scheme_dir"; then
        success=true
    fi

    # 标记完成
    if [ "$success" = true ]; then
        touch "$scheme_dir/.completed"
        echo "  ${GREEN}✓ 完成${NC}"
        return 0
    else
        echo "  ${RED}✗ 失败${NC}"
        return 1
    fi
}

# 主逻辑
main() {
    # 获取 schemes 列表
    schemes_json=$(get_schemes)

    if [ -z "$schemes_json" ]; then
        echo "${RED}错误: 无法获取 schemes 列表${NC}"
        exit 1
    fi

    # 检查是否安装了 jq
    if command -v jq &> /dev/null; then
        echo "使用 jq 解析 JSON"
        USE_JQ=1
        total=$(echo "$schemes_json" | jq '.schemes | length')
    else
        echo "未安装 jq，使用基本解析"
        USE_JQ=0
        total=$(echo "$schemes_json" | grep -o '"id"' | wc -l)
    fi

    echo "发现 $total 个 schemes"
    echo ""

    # 统计
    success_count=0
    failed_count=0
    skipped_count=0

    # 处理每个 scheme
    if [ "$USE_JQ" -eq 1 ]; then
        # 使用 jq 解析
        for i in $(seq 0 $((total - 1))); do
            scheme=$(echo "$schemes_json" | jq -r ".schemes[$i]")
            scheme_id=$(echo "$scheme" | jq -r '.id // .scheme_id')
            scheme_name=$(echo "$scheme" | jq -r '.name // "unknown"')
            species=$(echo "$scheme" | jq -r '.species // .organism // "unknown"')

            if process_scheme "$scheme_id" "$scheme_name" "$species" $((i + 1)) "$total"; then
                ((success_count++))
            else
                ((failed_count++))
            fi

            sleep "$DELAY"
        done
    else
        # 基本解析 - 仅支持简单的 scheme ID 提取
        echo "警告: 未安装 jq，仅支持基本功能"
        echo "建议安装 jq: apt-get install jq 或 brew install jq"

        # 尝试提取 ID 列表
        scheme_ids=$(echo "$schemes_json" | grep -o '"id":[0-9]*' | cut -d: -f2)

        index=1
        for scheme_id in $scheme_ids; do
            echo "[$index] Processing scheme ID: $scheme_id"

            scheme_dir="$OUTPUT_DIR/scheme_$scheme_id"
            mkdir -p "$scheme_dir"

            # 直接尝试下载
            download_scheme_detail "$scheme_id" "$scheme_dir"
            download_alleles "$scheme_id" "$scheme_dir"
            download_profiles "$scheme_id" "$scheme_dir"

            ((index++))
            sleep "$DELAY"
        done
    fi

    # 报告
    echo ""
    echo "========================================"
    echo "下载报告"
    echo "========================================"
    echo "总计: $total"
    echo "成功: $success_count"
    echo "失败: $failed_count"
    echo "跳过: $skipped_count"
    echo "输出目录: $(pwd)/$OUTPUT_DIR"
}

# 检查参数
if [ "$1" = "--test" ]; then
    echo "测试模式 - 仅获取 schemes 列表"
    get_schemes | head -c 2000
    echo "..."
    exit 0
fi

# 运行
main
