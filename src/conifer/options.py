"""Options for the CLI"""

import logging
from typing import Callable

import click
from click.decorators import FC

logger = logging.getLogger(__name__)


def verbose_option(expose_value: bool = False) -> Callable[[FC], FC]:
    """Add verbose option with callback"""

    def verbose_callback(
        ctx: click.core.Context, param: click.core.Option, value: int
    ) -> int:
        """Verbose callback"""
        log_level = max(3 - value, 0) * 10
        logging.basicConfig(
            level=log_level,
            format="%(levelname)s [%(name)s:%(funcName)s]: %(message)s",
        )
        return log_level

    verbose = click.option(
        "--verbose",
        "-v",
        help="Set the verbosity level",
        count=True,
        callback=verbose_callback,
        expose_value=expose_value,
        is_eager=True,
    )
    return verbose
