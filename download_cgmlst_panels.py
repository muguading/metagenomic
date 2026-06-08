#!/usr/bin/env python3
"""
cgMLST.org 全 Panel 下载脚本
支持下载所有物种的 cgMLST/wgMLST schemes
"""

import json
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Dict, List, Optional


BASE_URL = "https://www.cgmlst.org/ncs/api"
OUTPUT_DIR = Path("cgmlst_panels")

# 请求头
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "Accept": "application/json",
}


def fetch_json(url: str, max_retries: int = 3) -> Optional[dict]:
    """获取 JSON 数据，带重试机制"""
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=30) as response:
                return json.loads(response.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            if e.code == 429:  # Too Many Requests
                wait_time = (attempt + 1) * 5
                print(f"  请求过多，等待 {wait_time} 秒...")
                time.sleep(wait_time)
            else:
                print(f"  HTTP错误 {e.code}: {e.reason}")
                return None
        except Exception as e:
            print(f"  请求失败 (尝试 {attempt+1}/{max_retries}): {e}")
            time.sleep(2 ** attempt)
    return None


def fetch_text(url: str, max_retries: int = 3) -> Optional[str]:
    """获取文本数据，带重试机制"""
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=30) as response:
                return response.read().decode('utf-8')
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait_time = (attempt + 1) * 5
                print(f"  请求过多，等待 {wait_time} 秒...")
                time.sleep(wait_time)
            else:
                print(f"  HTTP错误 {e.code}: {e.reason}")
                return None
        except Exception as e:
            print(f"  请求失败 (尝试 {attempt+1}/{max_retries}): {e}")
            time.sleep(2 ** attempt)
    return None


def get_schemes_list() -> List[Dict]:
    """获取所有可用的 schemes 列表"""
    url = f"{BASE_URL}/schemes/"
    print(f"正在获取 schemes 列表: {url}")
    data = fetch_json(url)
    if data and "schemes" in data:
        return data["schemes"]
    return []


def get_scheme_details(scheme_id: int) -> Optional[Dict]:
    """获取单个 scheme 的详细信息"""
    url = f"{BASE_URL}/scheme/{scheme_id}"
    return fetch_json(url)


def download_alleles(scheme_id: int, scheme_name: str, output_dir: Path) -> bool:
    """下载所有 allele 序列 (FASTA格式)"""
    # 尝试多种可能的端点格式
    endpoints = [
        f"{BASE_URL}/scheme/{scheme_id}/alleles/fasta",
        f"{BASE_URL}/scheme/{scheme_id}/fasta",
        f"{BASE_URL}/alleles/{scheme_id}/fasta",
    ]

    for url in endpoints:
        print(f"  尝试下载 alleles: {url}")
        content = fetch_text(url)
        if content and len(content) > 100:
            output_file = output_dir / f"{scheme_name}_alleles.fasta"
            output_file.write_text(content, encoding='utf-8')
            print(f"  ✓ 已保存: {output_file}")
            return True
        time.sleep(1)

    return False


def download_profiles(scheme_id: int, scheme_name: str, output_dir: Path) -> bool:
    """下载 ST profiles (TSV格式)"""
    endpoints = [
        f"{BASE_URL}/scheme/{scheme_id}/profiles",
        f"{BASE_URL}/scheme/{scheme_id}/profiles.tsv",
        f"{BASE_URL}/profiles/{scheme_id}",
    ]

    for url in endpoints:
        print(f"  尝试下载 profiles: {url}")
        content = fetch_text(url)
        if content and len(content) > 50:
            output_file = output_dir / f"{scheme_name}_profiles.tsv"
            output_file.write_text(content, encoding='utf-8')
            print(f"  ✓ 已保存: {output_file}")
            return True
        time.sleep(1)

    return False


def download_loci(scheme_id: int, scheme_name: str, output_dir: Path) -> bool:
    """下载 loci/genes 信息"""
    url = f"{BASE_URL}/scheme/{scheme_id}"
    data = fetch_json(url)

    if data:
        output_file = output_dir / f"{scheme_name}_loci.json"
        output_file.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding='utf-8')
        print(f"  ✓ 已保存: {output_file}")
        return True
    return False


def safe_filename(name: str) -> str:
    """将名称转换为安全的文件名"""
    # 替换不安全字符
    unsafe_chars = '<>:"/\\|?* '
    for char in unsafe_chars:
        name = name.replace(char, '_')
    # 去除重复下划线
    while '__' in name:
        name = name.replace('__', '_')
    return name.strip('_')


def main():
    """主函数"""
    print("=" * 60)
    print("cgMLST.org 全 Panel 下载工具")
    print("=" * 60)

    # 创建输出目录
    OUTPUT_DIR.mkdir(exist_ok=True)
    print(f"输出目录: {OUTPUT_DIR.absolute()}")
    print()

    # 获取 schemes 列表
    schemes = get_schemes_list()
    if not schemes:
        print("错误: 无法获取 schemes 列表，请检查网络连接")
        sys.exit(1)

    print(f"发现 {len(schemes)} 个 schemes")
    print()

    # 下载统计
    stats = {
        "total": len(schemes),
        "success": 0,
        "failed": 0,
        "alleles_downloaded": 0,
        "profiles_downloaded": 0,
    }

    # 遍历下载每个 scheme
    for i, scheme in enumerate(schemes, 1):
        scheme_id = scheme.get("id", scheme.get("scheme_id", 0))
        scheme_name = scheme.get("name", scheme.get("scheme", f"scheme_{scheme_id}"))
        species = scheme.get("species", scheme.get("organism", "Unknown"))

        safe_name = safe_filename(f"{scheme_id}_{species}_{scheme_name}")

        print(f"[{i}/{len(schemes)}] 处理: {species} - {scheme_name} (ID: {scheme_id})")

        # 创建子目录
        scheme_dir = OUTPUT_DIR / safe_name
        scheme_dir.mkdir(exist_ok=True)

        # 保存 scheme 元数据
        metadata_file = scheme_dir / "metadata.json"
        metadata_file.write_text(json.dumps(scheme, indent=2, ensure_ascii=False), encoding='utf-8')

        # 下载各种数据
        success = False

        # 1. 下载 loci 信息
        if download_loci(scheme_id, "loci", scheme_dir):
            success = True

        # 2. 下载 alleles
        if download_alleles(scheme_id, "loci", scheme_dir):
            stats["alleles_downloaded"] += 1
            success = True
        else:
            print(f"  ⚠ 无法下载 alleles")

        # 3. 下载 profiles
        if download_profiles(scheme_id, "profiles", scheme_dir):
            stats["profiles_downloaded"] += 1
            success = True
        else:
            print(f"  ⚠ 无法下载 profiles")

        if success:
            stats["success"] += 1
        else:
            stats["failed"] += 1

        # 礼貌性延迟，避免服务器压力过大
        time.sleep(2)
        print()

    # 生成下载报告
    print("=" * 60)
    print("下载完成!")
    print("=" * 60)
    print(f"总计 schemes: {stats['total']}")
    print(f"成功下载: {stats['success']}")
    print(f"下载失败: {stats['failed']}")
    print(f"Alleles 下载: {stats['alleles_downloaded']}")
    print(f"Profiles 下载: {stats['profiles_downloaded']}")
    print(f"输出目录: {OUTPUT_DIR.absolute()}")

    # 保存下载日志
    log_file = OUTPUT_DIR / "download_log.json"
    log_file.write_text(json.dumps(stats, indent=2), encoding='utf-8')
    print(f"日志文件: {log_file}")


if __name__ == "__main__":
    main()
