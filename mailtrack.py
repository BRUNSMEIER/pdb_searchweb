import imaplib
import email
from email.header import decode_header
import os
import requests

# 你的邮箱和密码
EMAIL = "schonnbiukmit1@gmail.com"
PASSWORD = ""  # input ur password
IMAP_SERVER = "imap.gmail.com"

# 定义保存路径
output_dir = "downloaded_csvs"
os.makedirs(output_dir, exist_ok=True)

def clean(text):
    """清理字符串以生成有效文件名"""
    return "".join(c if c.isalnum() else "_" for c in text)

def download_csv(url, save_path):
    """下载 CSV 文件"""
    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, "wb") as file:
            file.write(response.content)
        print(f"Downloaded: {save_path}")
    else:
        print(f"Failed to download {url}. Status code: {response.status_code}")

def fetch_csv_links():
    """从邮箱中获取 CSV 链接"""
    csv_links = []

    # 连接到 Gmail 的 IMAP 服务器
    mail = imaplib.IMAP4_SSL(IMAP_SERVER)
    try:

        print(f"Attempting to log in with email: {EMAIL}")
        mail.login(EMAIL, PASSWORD)
        mail.select("inbox")
    except imaplib.IMAP4.error as e:
        print(f"IMAP error: {e}")


    # 搜索包含特定标题的邮件
    status, messages = mail.search(None, '(SUBJECT "predictions for nsite job id are ready")')

    # 获取邮件 ID
    for num in messages[0].split():
        # 获取邮件内容
        status, msg_data = mail.fetch(num, "(RFC822)")
        for response_part in msg_data:
            if isinstance(response_part, tuple):
                # 解析邮件内容
                msg = email.message_from_bytes(response_part[1])
                if msg.is_multipart():
                    for part in msg.walk():
                        if part.get_content_type() == "text/plain":  # 只解析纯文本部分
                            body = part.get_payload(decode=True).decode()
                            # 查找 CSV 链接
                            for line in body.splitlines():
                                if "The CSV file can be found here:" in line:
                                    csv_link = line.split("The CSV file can be found here: ")[-1]
                                    csv_links.append(csv_link)
                else:
                    body = msg.get_payload(decode=True).decode()
                    if "The CSV file can be found here:" in body:
                        csv_link = body.split("The CSV file can be found here: ")[-1]
                        csv_links.append(csv_link)

    # 退出邮箱
    mail.logout()
    return csv_links

def main():
    # 获取所有 CSV 链接
    csv_links = fetch_csv_links()
    print(f"Found {len(csv_links)} CSV links.")

    # 下载每个 CSV 文件
    for index, csv_link in enumerate(csv_links):
        filename = f"result_{index + 1}.csv"
        save_path = os.path.join(output_dir, filename)
        download_csv(csv_link, save_path)

if __name__ == "__main__":
    main()
