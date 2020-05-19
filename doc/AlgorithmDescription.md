<div style="font-size:3em; text-align:center;"> Algorithm Description </div>

# Abstract
This document describes technical aspects of coding tools included in
the associated codec. This document is not a specification of the associated
codec. Instead, it summarizes the highlighted features of coding tools for new
developers. This document should be updated when significant new normative
changes have been integrated into the associated codec.

# Abbreviations

CDEF: Constrained directional enhancement\
CfL: Chroma from Luma\
IntraBC: Intra block copy\
LCU: Largest coding unit\
OBMC: Overlapped Block Motion Compensation

# Algorithm Description

## Block Partitioning

### Coding block partition

The largest coding block unit (LCU) applied in this codec is 128×128. In
addition to no split mode `PARTITION_NONE`, the partition tree supports 9
different partitioning patterns, as shown in below figure.

<figure class="image"> <center><img src="img\partition_codingblock.svg"
alt="Partition" width="400" /> <figcaption>Figure 1: Supported coding block
partitions</figcaption> </figure>

According to the number of sub-partitions, the 9 partition modes are summarized
as follows: 1. Four partitions: `PARTITION_SPLIT`, `PARTITION_VERT_4`,
`PARTITION_HORZ_4` 2. Three partitions (T-Shape): `PARTITION_HORZ_A`,
`PARTITION_HORZ_B`, `PARTITION_VERT_A`, `PARTITION_HORZ_B` 3. Two partitions:
`PARTITION_HORZ`, `PARTITION_VERT`

Among all the 9 partitioning patterns, only `PARTITION_SPLIT` mode supports
recursive partitioning, i.e., sub-partitions can be further split, other
partitioning modes cannot further split. Particularly, for 8x8 and 128x128,
`PARTITION_VERT_4`, `PARTITION_HORZ_4` are not used, and for 8x8, T-Shape
partitions are not used either.

### Transform block partition

For inter prediction blocks, the coding block can be further partitioned into
multiple transform units. The transform unit partitioning can be done in a
recursive manner with the partitioning depth up to 2 levels. The transform
partitioning supports 1:1 (square), 1:2/2:1, and 1:4/4:1 transform unit sizes
ranging from 4×4 to 64×64. If the coding block is smaller than or equal to
64x64, the transform block partitioning can only apply to luma component, for
chroma blocks, the transform block size is identical to the coding block size.
Otherwise, if the coding block width or height is greater than 64, then both the
luma and chroma coding blocks will implicitly split into multiples of min(W,
64)x min(H, 64) and min(W, 32)x min(H, 32) transform blocks, respectively.

## Intra Prediction

### Directional intra prediction modes

Directional intra prediction modes are applied in intra prediction, which models
local textures using a given direction pattern. Directional intra prediction
modes are represented by nominal modes and angle delta. The nominal modes are
the same set of intra prediction angle used in VP9, which includes 8 angles. The
index value of angle delta is ranging from -3 ~ +3, and zero delta angle
indicates a nominal mode. The prediction angle is represented by a nominal intra
angle plus an angle delta. In total, there are 56 directional intra prediction
modes, as shown in the following figure. In the below figure, solid arrows
indicate directional intra prediction modes and dotted arrows represent non-zero
angle delta.

<figure class="image"> <center><img src="img\intra_directional.svg"
alt="Directional intra" width="300" /> <figcaption>Figure 2: Directional intra
prediction modes</figcaption> </figure>

The nominal mode index and angle delta index is signalled separately, and
nominal mode index is signalled before the associated angle delta index. It is
noted that for small block sizes, where the coding gain from extending intra
prediction angles may saturate, only the nominal modes are used and angle delta
index is not signalled.

### Non-directional intra prediction modes

In addition to directional intra prediction modes, four non-directional intra
modes which simulate smooth textures are also included. The four non-directional
intra modes include `SMOOTH_V`, `SMOOTH_H`, `SMOOTH` and `PAETH predictor`.

In `SMOOTH V`, `SMOOTH H` and `SMOOTH modes`, the prediction values are
generated using quadratic interpolation along vertical, horizontal directions,
or the average thereof. The samples used in the quadratic interpolation include
reconstructed samples from the top and left neighboring blocks and samples from
the right and bottom boundaries which are approximated by top reconstructed
samples and the left reconstructed samples.

In `PAETH predictor` mode, the prediction for each sample is assigned as one
from the top (T), left (L) and top-left (TL) reference samples, which has the
value closest to the Paeth predictor value, i.e., T + L -TL. The samples used in
`PAETH predictor` are illustrated in below figure.

<figure class="image"> <center><img src="img\intra_paeth.svg" alt="Directional
intra" width="300" /> <figcaption>Figure 3: Paeth predictor</figcaption>
</figure>

### Recursive filtering modes

Five filtering intra modes are defined, and each mode specify a set of eight
7-tap filters. Given the selected filtering mode index (0~4), the current block
is divided into 4x2 sub-blocks. For one 4×2 sub-block, each sample is predicted
by 7-tap interpolation using the 7 top and left neighboring samples as inputs.
Different filters are applied for samples located at different coordinates
within a 4×2 sub-block. The prediction process can be done recursively in unit
4x2 sub-block, which means that prediction samples generated for one 4x2
prediction block can be used to predict another 4x2 sub-block.

<figure class="image"> <center><img src="img\intra_recursive.svg"
alt="Directional intra" width="300" /> <figcaption>Figure 4: Recursive filtering
modes</figcaption> </figure>

### Chroma from Luma mode

Chroma from Luma (CfL) is a chroma intra prediction mode, which models chroma
samples as a linear function of co-located reconstructed luma samples. To align
the resolution between luma and chroma samples for different chroma sampling
format, e.g., 4:2:0 and 4:2:2, reconstructed luma pixels may need to be
sub-sampled before being used in CfL mode. In addition, the DC component is
removed to form the AC contribution. In CfL mode, the model parameters which
specify the linear function between two color components are optimized by
encoder signalled in the bitstream.

<figure class="image"> <center><img src="img\intra_cfl.svg"
alt="Directional intra" width="700" /> <figcaption>Figure 5: CfL
prediction</figcaption> </figure>

## Inter Prediction

### Motion vector prediction

Motion vectors are predicted by neighboring blocks which can be either spatial
neighboring blocks, or temporal neighboring blocks located in a reference frame.
A set of MV predictors will be identified by checking all these blocks and
utilized to encode the motion vector information.

**Spatial motion vector prediction**

There are two sets of spatial neighboring blocks that can be utilized for
finding spatial MV predictors, including the adjacent spatial neighbors which
are direct top and left neighbors of the current block, and second outer spatial
neighbors which are close but not directly adjacent to the current block. The
two sets of spatial neighboring blocks are illustrated in an example shown in
Figure 5.<figure class="image"> <center><img src="img\inter_spatial_mvp.svg"
alt="Directional intra" width="350" /><figcaption>Figure 6: Motion field
estimation by linear projection</figcaption></figure> For each set of spatial
neighbors, the top row will be checked from left to right and then the left
column will be checked from top to down. For the adjacent spatial neighbors, an
additional top-right block will be also checked after checking the left column
neighboring blocks. For the non-adjacent spatial neighbors, the top-left block
located at (-1, -1) position will be checked first, then the top row and left
column in a similar manner as the adjacent neighbors. The adjacent neighbors
will be checked first, then the temporal MV predictor that will be described in
the next subsection will be checked second, after that, the non-adjacent spatial
neighboring blocks will be checked.

For compound prediction which utilizes a pair of reference frames, the
non-adjacent spatial neighbors are not used for deriving the MV predictor.

**Temporal motion vector prediction**

In addition to spatial neighboring blocks, MV predictor can be also derived
using co-located blocks of reference pictures, namely temporal MV predictor. To
generate temporal MV predictor, the MVs of reference frames are first stored
together with reference indices associated with the reference frame. Then for
each 8x8 block of the current frame, the MVs of a reference frame which pass the
8x8 block are identified and stored together with the reference frame index in a
temporal MV buffer. In an example shown in Figure 5, the MV of reference frame 1
(R1) pointing from R1 to a reference frame of R1 is identified, i.e., MVref,
which passes a 8x8 block (shaded in blue dots) of current frame. Then this
MVref is stored in the temporal MV buffer associated with this 8x8 block.
<figure class="image"> <center><img src="img\inter_motion_field.svg"
alt="Directional intra" width="800" /><figcaption>Figure 7: Motion field
estimation by linear projection</figcaption></figure> Finally, given a couple of
pre-defined block coordinates, the associated MVs stored in the temporal MV
buffer are identified and projected accordingly to derive a temporal MV
predictor which points from the current block to its reference frame, e.g., MV0
in Figure 5. In Figure 6, the pre-defined block positions for deriving temporal
MV predictors of a 16x16 block are shown and up to 7 blocks will be checked to
find valid temporal MV predictors.<figure class="image"> <center><img
src="img\inter_tmvp_positions.svg" alt="Directional intra" width="300"
/><figcaption>Figure 8: Block positions for deriving temporal MV
predictors</figcaption></figure> The temporal MV predictors are checked after
the nearest spatial MV predictors but before the non-adjacent spatial MV
predictors.

All the spatial and temporal MV candidates will be put together in a pool, with
each predictor associated with a weighting determined during the scanning of the
spatial and temporal neighboring blocks. Based on the associated weightings, the
candidates are sorted and ranked, and up to four candidates will be used as a
list MV predictor list.

### Motion vector coding

### Interpolation filter for motion compensation

<mark>[Ed.: to be added]</mark>

### Warped motion compensation

**Global warped motion**

<mark>[Ed.: to be added]</mark>

**Local warped motion**

<mark>[Ed.: to be added]</mark>

### Overlapped block motion compensation

<mark>[Ed.: to be added]</mark>

### Reference frames

<mark>[Ed.: to be added]</mark>

### Compound Prediction

<mark>[Ed.: to be added]</mark>

**Compound wedge prediction**

<mark>[Ed.: to be added]</mark>

**Difference-modulated masked prediction**

<mark>[Ed.: to be added]</mark>

**Frame distance-based compound prediction**

<mark>[Ed.: to be added]</mark>

**Compound inter-intra prediction**

<mark>[Ed.: to be added]</mark>

## Transform

The separable 2D transform process is applied on prediction residuals. For the
forward transform, a 1-D vertical transform is performed first on each column of
the input residual block, then a horizontal transform is performed on each row
of the vertical transform output. For the backward transform, a 1-D horizontal
transform is performed first on each row of the input de-quantized coefficient
block, then a vertical transform is performed on each column of the horizontal
transform output. The primary 1-D transforms include four different types of
transform: a) 4-point, 8-point, 16-point, 32-point, 64-point DCT-2; b) 4-point,
8-point, 16-point asymmetric DST’s (DST-4, DST-7) and c) their flipped
versions; d) 4-point, 8-point, 16-point, 32-point identity transforms. When
transform size is 4-point, ADST refers to DST-7, otherwise, when transform size
is greater than 4-point, ADST refers to DST-4.


<center>

Table 1: Transform basis functions (DCT-2, DST-4 and DST-7 for N-point input.

| Transform Type          |      Basis function T<sub>_i_</sub>(_j_), _i_, _j_=0, 1,…, N−1 |
|-------------------------|:--------------------------------------------------------------:|
| DCT-2                   |  T<sub>_i_;</sub>(_j_) = &omega;<sub>0</sub> &sdot; _`sqrt`_(2/N) + _`cos`_( &pi;&sdot;_i_&sdot;(2&sdot;_j_ +1) ) / (2N) )<br>where &omega;<sub>0</sub> = (_i_==0) ? _`sqrt`_(N/2) : 1
| DST-4                   |    T<sub>_i_;</sub>(_j_) = _`sqrt`_(2/N) + _`sin`_( &pi;&sdot;(2&sdot;_i_+1)&sdot;(2&sdot;_j_+1) ) / (4&sdot;N) )                 |
| DST-7                   | T<sub>_i_;</sub>(_j_) = _`sqrt`_(4/(2&sdot;N+1)) + _`sin`_( &pi;&sdot;_i_&sdot;(_j_+1) ) / (2&sdot;N+1) ) |
| Identity (IDT)          | T<sub>_i_;</sub>(_j_) = (_i_==_j_) ? 1 : 0 |
|||
</center>

For luma component, each transform block can select one pair of horizontal and
vertical transform combination given a pre-defined set of transform type
candidates, and the selection is explicitly signalled into the bitstream. When
the maximum of transform block width and height is greater than or equal to 32,
the set of transform type candidates depend on the prediction mode, as described
in Table 2. Otherwise, when the maximum of transform block width and height is
smaller than 32, the set of transform type candidates depend on the prediction
mode, as described in Table 3.

<center>

Table 2: Transform type candidates for luma component when max(width, height)
is greater than or equal to 32.

| Max(width, height)      |      Intra                 |      Inter                    |
|:-----------------------:|:--------------------------:|:-----------------------------:|
| 32                      |  _DCTOnly_                 |      _DCTOnly_, _IDTX_        |
| 64                      |  _DCTOnly_                 |      _DCTOnly_                |
</center>

<center>

Table 3: Transform type candidates for luma component when max(width, height)
is smaller than 32.

| Min(width, height)      |      Intra                 |      Inter                    |
|:-----------------------:|:--------------------------:|:-----------------------------:|
| 4                       |  _DTT4_, _IDTX_, _1DDCT_   |      _ALL16_                  |
| 8                       |  _DTT4_, _IDTX_, _1DDCT_   |      _ALL16_                  |
| 16                      |  _DTT4_, _IDTX_            |      _DTT9, _IDTX_, _1DDCT_   |
</center>

The set of transform type candidates (namely transform set) is defined in
Table 3.

<center>

Table 3: Definition of transform set.

| Transform set | Vertical transform | Horizontal transform |
|---------------|:------------------:|:--------------------:|
| _DCTOnly_       |  DCT     |      DCT   |
| _IDTX_          |  IDT     |      IDT   |
| _1DDCT_         |  DCT <br> IDT              | IDT <br> DCT |
| _DTT4_          |  ADST <br> ADST <br> DCT <br> DCT              | ADST <br> DCT <br> ADST <br>  DCT                 |
| _DTT9_          |  DCT <br> DCT <br> DCT <br> ADST <br> ADST <br> ADST <br> Flipped ADST <br> Flipped ADST <br> Flipped ADST              | DCT <br> ADST <br> Flipped ADST <br> DCT <br> ADST <br> Flipped ADST <br> DCT <br> ADST <br> Flipped ADST <br>                 |
| _ALL16_          |  DCT <br> DCT <br> DCT <br> DCT <br> ADST <br> ADST <br>ADST <br> ADST <br> Flipped ADST <br> Flipped ADST <br> Flipped ADST <br> Flipped ADST <br> IDT <br> IDT <br> IDT <br> IDT <br>          | DCT <br> ADST <br> Flipped ADST <br> IDT <br> DCT <br> ADST <br> Flipped ADST <br> IDT <br> DCT <br> ADST <br> Flipped ADST <br> IDT <br> DCT <br> ADST <br> Flipped ADST <br> IDT <br>                 |
||
</center>

For chroma component, the transform type selection is done in an implicit way.
For intra prediction residuals, the transform type is selected according to the
intra prediction mode, as specified in Table 4. For inter prediction residuals,
the transform type is selected according to the transform type selection of the
co-located luma block. Therefore, for chroma component, there is no transform
type signalling in the bitstream.

<center>

Table 4: Transform type selection for chroma component intra prediction residuals.


| Intra Prediction Mode | Vertical transform | Horizontal transform |
|-----------------------|:------------------:|:--------------------:|
| DC_PRED               | DCT                | DCT                  |
| V_PRED                | ADST               | DCT                  |
| H_PRED                | DCT                | ADST                 |
| D45_PRED              | DCT                | DCT                  |
| D135_PRED             | ADST               | ADST                 |
| D113_PRED             | ADST               | DCT                  |
| D157_PRED             | DCT                | ADST                 |
| D203_PRED             | DCT                | ADST                 |
| D67_PRED              | ADST               | DCT                  |
| SMOOTH_PRED           | ADST               | ADST                 |
| SMOOTH_V_PRED         | ADST               | DCT                  |
| SMOOTH_H_PRED         | DCT                | ADST                 |
| PAETH_PRED            | ADST               | ADST                 |
||
</center>

The computational cost of large size (e.g., 64-point) transforms is further
reduced by zeroing out all the coefficients except the following two cases:

1. The top-left 32×32 quadrant for 64×64/64×32/32×64 DCT_DCT hybrid transforms
2. The top-left 16×16 quadrant for 64×16/16×64 DCT_DCT hybrid transforms.

Both the DCT-2 and ADST (DST-4, DST-7) are implemented using butterfly structure
[1], which included multiple stages of butterfly operations. Each butterfly
operations can be calculated in parallel and different stages are cascaded in a
sequential order.

## Quantization
<mark>[Ed.: to be added]</mark>

## Entropy Coding

**Entropy coding engine**

<mark>[Ed.: to be added]</mark>

**Coefficient coding**

For each transform unit, the coefficient coding starts with coding a skip sign,
which is followed by the signaling of primary transform kernel type and the
end-of-block (EOB) position in case the transform coding is not skipped. After
that, the coefficient values are coded in a multiple level map manner plus sign
values. The level maps are coded as three level planes, namely lower-level,
middle-level and higher-level planes, and the sign is coded as another separate
plane. The lower-level, middle-level and higher-level planes correspond to
correspond to different ranges of coefficient magnitudes. The lower level plane
corresponds to the range of 0–2, the middle level plane takes care of the
range of 3–14, and the higher-level plane covers the range of 15 and above.

The three level planes are coded as follows. After the EOB position is coded,
the lower-level and middle-level planes are coded together in backward scan
order, and the scan order refers to zig-zag scan applied on the entire transform
unit basis. Then the sign plane and higher-level plane are coded together in
forward scan order. After that, the remainder (coefficient level minus 14) is
entropy coded using Exp-Golomb code.

The context model applied to the lower level plane depends on the primary
transform directions, including: bi-directional, horizontal, and vertical, as
well as transform size, and up to five neighbor (in frequency domain)
coefficients are used to derive the context. The middle level plane uses a
similar context model, but the number of context neighbor coefficients is
reduced from 5 to 2. The higher-level plane is coded by Exp-Golomb code without
using context model. For the sign plane, except the DC sign that is coded using
the DC signs from its neighboring transform units, sign values of other
coefficients are coded directly without using context model.

## Loop filtering and post-processing

### Deblocking

<mark>[Ed.: to be added]</mark>

### Constrained directional enhancement

<mark>[Ed.: to be added]</mark>

**Direction Estimation**

<mark>[Ed.: to be added]</mark>

**Non-linear low-pass filter**

<mark>[Ed.: to be added]</mark>

### Loop Restoration filter

**Separable symmetric normalized Wiener filter**

<mark>[Ed.: to be added]</mark>

**Dual self-guided filter**

<mark>[Ed.: to be added]</mark>

### Frame super-resolution

<mark>[Ed.: to be added]</mark>

### Film grain synthesis

<mark>[Ed.: to be added]</mark>

## Screen content coding

### Intra block copy

<mark>[Ed.: to be added]</mark>

### Palette mode

<mark>[Ed.: to be added]</mark>

# References

[1] J. Han, Y. Xu and D. Mukherjee, "A butterfly structured design of the hybrid
transform coding scheme," 2013 Picture Coding Symposium (PCS), San Jose, CA,
2013, pp. 17-20.