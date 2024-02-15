from masq.utils.regions import Region
from masq.utils.regions import collapse_regions


def test_region_from_string() -> None:
    region = "chr1:100-200"
    result = Region.from_string(region)
    assert result.chrom == "chr1"
    assert result.start == 100
    assert result.stop == 200


def test_collapse_regions() -> None:
    regions = [
        Region("chr1", 100, 200),
        Region("chr1", 150, 250),
        Region("chr1", 300, 400),
    ]
    result = collapse_regions(regions)
    assert result == [Region("chr1", 100, 250), Region("chr1", 300, 400)]
