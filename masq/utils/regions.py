from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Region:
    chrom: str
    start: int
    stop: int

    @staticmethod
    def from_string(region: str) -> Region:
        chrom, start_stop = region.split(":")
        start, stop = map(int, start_stop.split("-"))
        assert start <= stop
        return Region(chrom, start, stop)


def collapse_regions(regions: list[Region]) -> list[Region]:
    regions.sort(key=lambda x: (x.chrom, x.start))
    collapsed = []
    current = regions[0]
    for region in regions[1:]:
        if region.chrom != current.chrom:
            collapsed.append(current)
            current = region
        elif region.start <= current.stop:
            current = Region(
                current.chrom, current.start,
                max(current.stop, region.stop))
        else:
            collapsed.append(current)
            current = region
    collapsed.append(current)
    return collapsed
