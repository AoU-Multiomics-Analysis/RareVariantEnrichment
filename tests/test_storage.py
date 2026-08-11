from pathlib import Path

import pytest

from rare_variant_enrichment.storage import MinimumDistanceStore


def test_store_reduces_within_but_not_across_annotation_classes(tmp_path: Path):
    with MinimumDistanceStore(tmp_path) as store:
        store.upsert("S1", "ENSG1.1", "AC=1", "consequence", "stop_gained", 50)
        store.upsert("S1", "ENSG1.1", "AC=1", "consequence", "missense_variant", 20)
        store.upsert("S1", "ENSG1.1", "AC=1", "consequence", "stop_gained", 10)

        assert list(store.iter_feature("ENSG1.1")) == [
            ("S1", "AC=1", "consequence", "missense_variant", 20),
            ("S1", "AC=1", "consequence", "stop_gained", 10),
        ]


@pytest.mark.parametrize("annotation_family, annotation_class", [("", "stop_gained"), ("consequence", "")])
def test_store_rejects_empty_annotation_key_fields(
    tmp_path: Path, annotation_family: str, annotation_class: str
):
    with MinimumDistanceStore(tmp_path) as store:
        with pytest.raises(ValueError, match="annotation family and class must be non-empty"):
            store.upsert(
                "S1", "ENSG1.1", "AC=1", annotation_family, annotation_class, 10
            )
