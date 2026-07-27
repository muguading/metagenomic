cd /data/HPCDC_pip/linux_minimal_portal

./.venv_web/bin/python - <<'PY'
import sqlite3
from werkzeug.security import generate_password_hash

username = "admin"
new_password = "admin123"

conn = sqlite3.connect("bac_analysis_portal.sqlite3")
cur = conn.cursor()
cur.execute(
    "UPDATE users SET password_hash=? WHERE username=?",
    (generate_password_hash(new_password), username),
)
conn.commit()
print("updated rows:", cur.rowcount)
conn.close()
PY