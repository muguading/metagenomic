from __future__ import annotations

from .application import create_app
from .knowledge_interpretation import _build_viral_serotype_knowledge_summary


if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5055, debug=app.config["PORTAL_MODE"] == "development")
