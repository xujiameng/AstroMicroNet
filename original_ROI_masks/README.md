# Original microdomain ROI masks

Thirteen astrocytes from five mice, with 3,236 unique microdomain masks mapped to 3,291 original records. These archived masks were copied without redraw or reassignment. The export was previously checked against the original ROIInfo.mask arrays; the present package preserves the MAT and NPZ bytes.

`cells.csv` identifies the cells and mice. `rois.csv` gives original record indices, areas and mask hashes. In each MAT/NPZ file, sparse rows are ROI masks and columns use MATLAB column-major pixel order; fields ending in `_1based` are one-based. `image_shape` must be read for each cell. `source_record_index_1based` and `mapback_1based` retain the original mapping. Branch territory masks describe original manualMask territories, not the excluded soma/process/endfoot masks.

The submitted procedure delineates the cell in the mean structural image, excludes soma, primary processes and endfeet, and uses local morphology and fluorescence to delineate remaining microdomains. These are qualitative historical rules; no new uniform area cutoff is imposed on the archived masks. Their actual areas and overlaps remain unchanged. The sampling-comparison ROIs and public H1R ROIs are separate datasets.

MATLAB: `S=load('20210805A07_1.mat'); imshow(reshape(full(S.roi_masks(1,:)),S.image_shape(1),S.image_shape(2)));`
