import pytest

from itertools import product

from utils.helpers import fetch_main_targets
from .helpers import run_dag

targets = fetch_main_targets()
print(targets)

remove_related = [True, False]


@pytest.mark.dryrun
@pytest.mark.parametrize("target,remove_related", product(targets, remove_related))
def test_saige_dag(target, remove_related):

    remove_related = "remove_related={}".format(remove_related)
    sample_sheet = "sample_sheet={}".format("tests/test_data/sample_sheet.txt")
    config = "--config {}".format(remove_related)

    exit_status = run_dag(target, "tests/test_data/saige_example/config.yaml", config=config)

    assert exit_status == 0

# because king not in bioconda yet
remove_related_integration = [False]

@pytest.mark.integration
@pytest.mark.parametrize("target,remove_related_integration", product(targets, remove_related_integration))
def test_saige_example(target, remove_related_integration):

    remove_related = "remove_related={}".format(remove_related_integration)

    exit_status = run_dag(target, "tests/test_data/saige_example/config.yaml", config=remove_related, dryrun=False,
                          extras="--use-conda")

    assert exit_status == 0
