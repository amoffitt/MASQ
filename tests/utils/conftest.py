import os
from typing import Callable

import pytest


@pytest.fixture(scope='session')
def fixtures_filename() -> Callable[[str], str]:

    def build_filename(filename: str) -> str:
        return os.path.join(
            os.path.dirname(__file__),
            "fixtures",
            filename)

    return build_filename
