from bs4 import BeautifulSoup
import os


def parse_results(result_file):
    with open(result_file, "r") as f:
        soup = BeautifulSoup(f, "html.parser")

    # 假设结果在某个表格中
    table = soup.find("table", {"id": "resultsTable"})
    rows = table.find_all("tr")[1:]  # 跳过表头
    results = []
    for row in rows:
        cols = row.find_all("td")
        data = [col.text.strip() for col in cols]
        results.append(data)
    return results


# 遍历结果文件
result_files = [f for f in os.listdir() if f.endswith("_result.html")]
for result_file in result_files:
    data = parse_results(result_file)
    print(f"Parsed data for {result_file}: {data}")
