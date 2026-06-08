#!/usr/bin/env python3
"""
测试 cgMLST.org API 连接性和端点结构
"""

import json
import urllib.request
import urllib.error


def test_endpoint(base_url: str) -> bool:
    """测试单个端点"""
    print(f"\n测试端点: {base_url}")
    print("-" * 50)

    schemes_url = f"{base_url}/schemes/"
    print(f"尝试: {schemes_url}")

    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
        "Accept": "application/json",
    }

    try:
        req = urllib.request.Request(schemes_url, headers=headers)
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read().decode('utf-8')
            content_type = resp.headers.get('Content-Type', '')
            print(f"  状态: {resp.status}")
            print(f"  内容类型: {content_type}")
            print(f"  数据长度: {len(data)} bytes")

            # 尝试解析 JSON
            try:
                json_data = json.loads(data)
                print(f"  JSON 结构: {type(json_data).__name__}")

                if isinstance(json_data, dict):
                    print(f"  顶层键: {list(json_data.keys())}")
                    if "schemes" in json_data:
                        schemes = json_data["schemes"]
                        print(f"  schemes 数量: {len(schemes)}")
                        if schemes:
                            print(f"  第一个 scheme 结构:")
                            first = schemes[0]
                            for key, val in first.items():
                                print(f"    {key}: {type(val).__name__} = {val}")

                            # 测试获取单个 scheme 详情
                            scheme_id = first.get("id", first.get("scheme_id"))
                            if scheme_id:
                                test_scheme_detail(base_url, scheme_id, headers)

                return True

            except json.JSONDecodeError as e:
                print(f"  JSON 解析失败: {e}")
                print(f"  原始数据前 500 字符: {data[:500]}")
                return False

    except urllib.error.HTTPError as e:
        print(f"  HTTP 错误: {e.code} - {e.reason}")
        return False
    except Exception as e:
        print(f"  错误: {e}")
        return False


def test_scheme_detail(base_url: str, scheme_id: int, headers: dict):
    """测试获取 scheme 详情"""
    print(f"\n  测试 scheme 详情 (ID: {scheme_id})")
    detail_url = f"{base_url}/scheme/{scheme_id}"
    print(f"  尝试: {detail_url}")

    try:
        req = urllib.request.Request(detail_url, headers=headers)
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            print(f"    成功! 键: {list(data.keys())}")

            # 测试下载链接
            test_downloads(base_url, scheme_id, headers)

    except Exception as e:
        print(f"    失败: {e}")


def test_downloads(base_url: str, scheme_id: int, headers: dict):
    """测试数据下载链接"""
    print(f"\n  测试数据下载链接 (ID: {scheme_id})")

    endpoints = [
        ("alleles/fasta", f"{base_url}/scheme/{scheme_id}/alleles/fasta"),
        ("fasta", f"{base_url}/scheme/{scheme_id}/fasta"),
        ("profiles", f"{base_url}/scheme/{scheme_id}/profiles"),
        ("profiles.tsv", f"{base_url}/scheme/{scheme_id}/profiles.tsv"),
    ]

    for name, url in endpoints:
        print(f"    测试 {name}: ", end="")
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = resp.read()
                if len(data) > 100:
                    preview = data[:100].decode('utf-8', errors='ignore').replace('\n', ' ')
                    print(f"✓ ({len(data):,} bytes) - {preview}...")
                else:
                    print(f"✗ (数据太短: {len(data)} bytes)")
        except urllib.error.HTTPError as e:
            print(f"✗ HTTP {e.code}")
        except Exception as e:
            print(f"✗ {e}")


def main():
    """主函数"""
    print("=" * 60)
    print("cgMLST.org API 测试工具")
    print("=" * 60)

    # 测试多个可能的端点
    endpoints = [
        "https://www.cgmlst.org/ncs/api",
        "https://cgmlst.org/ncs/api",
        "https://www.cgmlst.org/api",
        "https://cgmlst.org/api",
    ]

    for endpoint in endpoints:
        if test_endpoint(endpoint):
            print(f"\n✓ 可用端点: {endpoint}")
            break
    else:
        print("\n✗ 所有端点均不可用")

        print("\n可能的原因:")
        print("  1. 网站需要浏览器级别的 JavaScript 支持")
        print("  2. 网站有反爬虫机制")
        print("  3. 网络连接问题")
        print("  4. API 端点已变更")
        print("\n建议:")
        print("  1. 使用浏览器开发者工具(F12)查看实际 API 请求")
        print("  2. 检查网站是否有 robots.txt 限制")
        print("  3. 尝试使用浏览器插件导出 HAR 文件分析")


if __name__ == "__main__":
    main()
