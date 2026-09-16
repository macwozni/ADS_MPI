"""Allow ``python -m ads_benchmark`` to invoke the planner."""

from .cli import main


raise SystemExit(main())
