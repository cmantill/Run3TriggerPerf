Ë² 
§))
,
Abs
x"T
y"T"
Ttype:

2	
:
Add
x"T
y"T
z"T"
Ttype:
2	
x
Assign
ref"T

value"T

output_ref"T"	
Ttype"
validate_shapebool("
use_lockingbool(
s
	AssignSub
ref"T

value"T

output_ref"T" 
Ttype:
2	"
use_lockingbool( 
k
BatchMatMulV2
x"T
y"T
output"T"
Ttype:

2	"
adj_xbool( "
adj_ybool( 
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
N
Cast	
x"SrcT	
y"DstT"
SrcTtype"
DstTtype"
Truncatebool( 
h
ConcatV2
values"T*N
axis"Tidx
output"T"
Nint(0"	
Ttype"
Tidxtype0:
2	
8
Const
output"dtype"
valuetensor"
dtypetype

Conv2D

input"T
filter"T
output"T"
Ttype:
2"
strides	list(int)"
use_cudnn_on_gpubool(",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 "-
data_formatstringNHWC:
NHWCNCHW" 
	dilations	list(int)

B
Equal
x"T
y"T
z
"
Ttype:
2	

W

ExpandDims

input"T
dim"Tdim
output"T"	
Ttype"
Tdimtype0:
2	
^
Fill
dims"
index_type

value"T
output"T"	
Ttype"

index_typetype0:
2	
­
GatherV2
params"Tparams
indices"Tindices
axis"Taxis
output"Tparams"

batch_dimsint "
Tparamstype"
Tindicestype:
2	"
Taxistype:
2	
B
GreaterEqual
x"T
y"T
z
"
Ttype:
2	
.
Identity

input"T
output"T"	
Ttype
2
L2Loss
t"T
output"T"
Ttype:
2
\
	LeakyRelu
features"T
activations"T"
alphafloat%ÍÌL>"
Ttype0:
2
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	

Max

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
Ô
MaxPool

input"T
output"T"
Ttype0:
2	"
ksize	list(int)(0"
strides	list(int)(0""
paddingstring:
SAMEVALID":
data_formatstringNHWC:
NHWCNCHWNCHW_VECT_C

Mean

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
N
Merge
inputs"T*N
output"T
value_index"	
Ttype"
Nint(0
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(
=
Mul
x"T
y"T
z"T"
Ttype:
2	
.
Neg
x"T
y"T"
Ttype:

2	

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
X
PlaceholderWithDefault
input"dtype
output"dtype"
dtypetype"
shapeshape
~
RandomUniform

shape"T
output"dtype"
seedint "
seed2int "
dtypetype:
2"
Ttype:
2	
a
Range
start"Tidx
limit"Tidx
delta"Tidx
output"Tidx"
Tidxtype0:	
2	
>
RealDiv
x"T
y"T
z"T"
Ttype:
2	
\
	RefSwitch
data"T
pred

output_false"T
output_true"T"	
Ttype
E
Relu
features"T
activations"T"
Ttype:
2	
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
.
Rsqrt
x"T
y"T"
Ttype:

2
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
9
Softmax
logits"T
softmax"T"
Ttype:
2
1
Square
x"T
y"T"
Ttype:

2	
G
SquaredDifference
x"T
y"T
z"T"
Ttype:

2	
N
Squeeze

input"T
output"T"	
Ttype"
squeeze_dims	list(int)
 (
2
StopGradient

input"T
output"T"	
Ttype
ö
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
:
Sub
x"T
y"T
z"T"
Ttype:
2	

Sum

input"T
reduction_indices"Tidx
output"T"
	keep_dimsbool( " 
Ttype:
2	"
Tidxtype0:
2	
M
Switch	
data"T
pred

output_false"T
output_true"T"	
Ttype
c
Tile

input"T
	multiples"
Tmultiples
output"T"	
Ttype"

Tmultiplestype0:
2	
f
TopKV2

input"T
k
values"T
indices"
sortedbool("
Ttype:
2	
P
	Transpose
x"T
perm"Tperm
y"T"	
Ttype"
Tpermtype0:
2	
s

VariableV2
ref"dtype"
shapeshape"
dtypetype"
	containerstring "
shared_namestring "serve*1.14.12v1.14.0-65-g00fad90®
d
PlaceholderPlaceholder*
dtype0*"
_output_shapes
:*
shape:
V
Placeholder_1Placeholder*
_output_shapes
:*
shape:*
dtype0
N
Placeholder_2Placeholder*
dtype0
*
_output_shapes
: *
shape: 
Y
ExpandDims/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
r

ExpandDims
ExpandDimsPlaceholderExpandDims/dim*
T0*&
_output_shapes
:*

Tdim0
h
strided_slice/stackConst*!
valueB"            *
dtype0*
_output_shapes
:
j
strided_slice/stack_1Const*!
valueB"           *
dtype0*
_output_shapes
:
j
strided_slice/stack_2Const*!
valueB"         *
dtype0*
_output_shapes
:

strided_sliceStridedSlicePlaceholderstrided_slice/stackstrided_slice/stack_1strided_slice/stack_2*
shrink_axis_mask *
ellipsis_mask *

begin_mask*
new_axis_mask *
end_mask*"
_output_shapes
:*
T0*
Index0
^
SqueezeSqueezestrided_slice*
T0*
_output_shapes

:*
squeeze_dims
 
R
ExpandDims_1/dimConst*
value	B : *
dtype0*
_output_shapes
: 
n
ExpandDims_1
ExpandDimsSqueezeExpandDims_1/dim*"
_output_shapes
:*

Tdim0*
T0
j
strided_slice_1/stackConst*!
valueB"           *
dtype0*
_output_shapes
:
l
strided_slice_1/stack_1Const*!
valueB"           *
dtype0*
_output_shapes
:
l
strided_slice_1/stack_2Const*
_output_shapes
:*!
valueB"         *
dtype0

strided_slice_1StridedSliceExpandDims_1strided_slice_1/stackstrided_slice_1/stack_1strided_slice_1/stack_2*
Index0*
T0*
shrink_axis_mask*

begin_mask*
ellipsis_mask *
new_axis_mask *
end_mask*
_output_shapes

:
[
ExpandDims_2/dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
v
ExpandDims_2
ExpandDimsstrided_slice_1ExpandDims_2/dim*"
_output_shapes
:*

Tdim0*
T0
f
strided_slice_2/stackConst*
valueB"        *
dtype0*
_output_shapes
:
h
strided_slice_2/stack_1Const*
valueB"        *
dtype0*
_output_shapes
:
h
strided_slice_2/stack_2Const*
_output_shapes
:*
valueB"      *
dtype0

strided_slice_2StridedSliceExpandDims_1strided_slice_2/stackstrided_slice_2/stack_1strided_slice_2/stack_2*
Index0*
T0*
shrink_axis_mask *
ellipsis_mask *

begin_mask*
new_axis_mask *
end_mask*"
_output_shapes
:
L
Equal/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
R
EqualEqualExpandDims_2Equal/y*"
_output_shapes
:*
T0
d
ones_like/ShapeConst*!
valueB"         *
dtype0*
_output_shapes
:
T
ones_like/ConstConst*
_output_shapes
: *
valueB
 *  ?*
dtype0
r
	ones_likeFillones_like/Shapeones_like/Const*
T0*

index_type0*"
_output_shapes
:
Z
ShapeConst*
_output_shapes
:*!
valueB"         *
dtype0
J
ConstConst*
valueB
 *    *
dtype0*
_output_shapes
: 
Y
FillFillShapeConst*
T0*

index_type0*"
_output_shapes
:
U
SelectSelectEqual	ones_likeFill*
T0*"
_output_shapes
:
J
mul/xConst*
_output_shapes
: *
valueB
 *  zD*
dtype0
F
mulMulmul/xSelect*
T0*"
_output_shapes
:
c
transpose/permConst*!
valueB"          *
dtype0*
_output_shapes
:
e
	transpose	Transposemultranspose/perm*"
_output_shapes
:*
Tperm0*
T0
j
strided_slice_3/stackConst*
_output_shapes
:*!
valueB"            *
dtype0
l
strided_slice_3/stack_1Const*
_output_shapes
:*!
valueB"           *
dtype0
l
strided_slice_3/stack_2Const*!
valueB"         *
dtype0*
_output_shapes
:

strided_slice_3StridedSliceExpandDims_1strided_slice_3/stackstrided_slice_3/stack_1strided_slice_3/stack_2*
shrink_axis_mask *

begin_mask*
ellipsis_mask *
new_axis_mask *
end_mask*"
_output_shapes
:*
Index0*
T0
e
transpose_1/permConst*!
valueB"          *
dtype0*
_output_shapes
:
u
transpose_1	Transposestrided_slice_3transpose_1/perm*
T0*"
_output_shapes
:*
Tperm0
j
strided_slice_4/stackConst*!
valueB"           *
dtype0*
_output_shapes
:
l
strided_slice_4/stack_1Const*
_output_shapes
:*!
valueB"            *
dtype0
l
strided_slice_4/stack_2Const*!
valueB"         *
dtype0*
_output_shapes
:

strided_slice_4StridedSlicetranspose_1strided_slice_4/stackstrided_slice_4/stack_1strided_slice_4/stack_2*
ellipsis_mask *

begin_mask*
new_axis_mask *
end_mask*"
_output_shapes
:*
T0*
Index0*
shrink_axis_mask 
c
Tile/multiplesConst*!
valueB"         *
dtype0*
_output_shapes
:
l
TileTilestrided_slice_4Tile/multiples*"
_output_shapes
:*

Tmultiples0*
T0
e
transpose_2/permConst*!
valueB"          *
dtype0*
_output_shapes
:
j
transpose_2	TransposeTiletranspose_2/perm*
T0*"
_output_shapes
:*
Tperm0
J
subSubTiletranspose_2*"
_output_shapes
:*
T0
<
AbsAbssub*
T0*"
_output_shapes
:
>
Abs_1AbsAbs*
T0*"
_output_shapes
:
S
GreaterEqual/yConst*
valueB
 *ÛÉ@*
dtype0*
_output_shapes
: 
`
GreaterEqualGreaterEqualAbs_1GreaterEqual/y*"
_output_shapes
:*
T0
L
mul_1/xConst*
valueB
 *ÛIA*
dtype0*
_output_shapes
: 
G
mul_1Mulmul_1/xAbs*"
_output_shapes
:*
T0
L
sub_1/xConst*
valueB
 *æéB*
dtype0*
_output_shapes
: 
I
sub_1Subsub_1/xmul_1*
T0*"
_output_shapes
:
C
sub_2SubAbsAbs*
T0*"
_output_shapes
:
[
Select_1SelectGreaterEqualsub_1sub_2*
T0*"
_output_shapes
:
|
MatMulBatchMatMulV2strided_slice_3transpose_1*
adj_y( *
T0*"
_output_shapes
:*
adj_x( 
L
mul_2/xConst*
_output_shapes
: *
valueB
 *   À*
dtype0
J
mul_2Mulmul_2/xMatMul*
T0*"
_output_shapes
:
N
SquareSquarestrided_slice_3*"
_output_shapes
:*
T0
`
Sum/reduction_indicesConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
s
SumSumSquareSum/reduction_indices*"
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
e
transpose_3/permConst*
_output_shapes
:*!
valueB"          *
dtype0
i
transpose_3	TransposeSumtranspose_3/perm*"
_output_shapes
:*
Tperm0*
T0
C
addAddSummul_2*
T0*"
_output_shapes
:
K
add_1Addaddtranspose_3*"
_output_shapes
:*
T0
J
add_2Addadd_1Select_1*
T0*"
_output_shapes
:
E
add_3Addadd_2mul*"
_output_shapes
:*
T0
K
add_4Addadd_3	transpose*
T0*"
_output_shapes
:
>
NegNegadd_4*
T0*"
_output_shapes
:
J
TopKV2/kConst*
_output_shapes
: *
value	B :
*
dtype0
h
TopKV2TopKV2NegTopKV2/k*
sorted(*
T0*0
_output_shapes
:
:

^
	Squeeze_1SqueezePlaceholder*
T0*
_output_shapes

:*
squeeze_dims
 
R
ExpandDims_3/dimConst*
value	B : *
dtype0*
_output_shapes
: 
p
ExpandDims_3
ExpandDims	Squeeze_1ExpandDims_3/dim*

Tdim0*
T0*"
_output_shapes
:
[
ExpandDims_4/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
w
ExpandDims_4
ExpandDimsExpandDims_3ExpandDims_4/dim*&
_output_shapes
:*

Tdim0*
T0
á
Jlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/shapeConst*
_output_shapes
:*%
valueB"            *<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
dtype0
Ë
Hlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *Bãè¾*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ë
Hlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/maxConst*
_output_shapes
: *
valueB
 *Bãè>*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
dtype0
Ä
Rlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/shape*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
seed2 *
dtype0*&
_output_shapes
:*

seed 
Â
Hlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/subSubHlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/maxHlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
_output_shapes
: 
Ü
Hlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/mulMulRlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformHlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/sub*&
_output_shapes
:*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights
Î
Dlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniformAddHlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/mulHlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform/min*&
_output_shapes
:*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights
ú
)layerfilter0PL_newfea_conv_head_0/weights
VariableV2"/device:CPU:0*
shared_name *<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
	container *
shape:*
dtype0*&
_output_shapes
:
Ò
0layerfilter0PL_newfea_conv_head_0/weights/AssignAssign)layerfilter0PL_newfea_conv_head_0/weightsDlayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*&
_output_shapes
:*
use_locking(*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
validate_shape(
ã
.layerfilter0PL_newfea_conv_head_0/weights/readIdentity)layerfilter0PL_newfea_conv_head_0/weights"/device:CPU:0*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*&
_output_shapes
:

(layerfilter0PL_newfea_conv_head_0/L2LossL2Loss.layerfilter0PL_newfea_conv_head_0/weights/read*
_output_shapes
: *
T0
t
/layerfilter0PL_newfea_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter0PL_newfea_conv_head_0/weight_lossMul(layerfilter0PL_newfea_conv_head_0/L2Loss/layerfilter0PL_newfea_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 
©
(layerfilter0PL_newfea_conv_head_0/Conv2DConv2DExpandDims_4.layerfilter0PL_newfea_conv_head_0/weights/read*&
_output_shapes
:*
	dilations
*
T0*
data_formatNHWC*
strides
*
explicit_paddings
 *
use_cudnn_on_gpu(*
paddingVALID
w
*layerfilter0PL_newfea_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

)layerfilter0PL_newfea_conv_head_0/bn/beta
VariableV2*
_output_shapes
:*
	container *
shape:*
shared_name *
dtype0

0layerfilter0PL_newfea_conv_head_0/bn/beta/AssignAssign)layerfilter0PL_newfea_conv_head_0/bn/beta*layerfilter0PL_newfea_conv_head_0/bn/Const*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
È
.layerfilter0PL_newfea_conv_head_0/bn/beta/readIdentity)layerfilter0PL_newfea_conv_head_0/bn/beta*
_output_shapes
:*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/bn/beta
y
,layerfilter0PL_newfea_conv_head_0/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

*layerfilter0PL_newfea_conv_head_0/bn/gamma
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes
:*
	container 
¢
1layerfilter0PL_newfea_conv_head_0/bn/gamma/AssignAssign*layerfilter0PL_newfea_conv_head_0/bn/gamma,layerfilter0PL_newfea_conv_head_0/bn/Const_1*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_newfea_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
Ë
/layerfilter0PL_newfea_conv_head_0/bn/gamma/readIdentity*layerfilter0PL_newfea_conv_head_0/bn/gamma*=
_class3
1/loc:@layerfilter0PL_newfea_conv_head_0/bn/gamma*
_output_shapes
:*
T0

Clayerfilter0PL_newfea_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ö
1layerfilter0PL_newfea_conv_head_0/bn/moments/meanMean(layerfilter0PL_newfea_conv_head_0/Conv2DClayerfilter0PL_newfea_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*

Tidx0*
T0*&
_output_shapes
:
­
9layerfilter0PL_newfea_conv_head_0/bn/moments/StopGradientStopGradient1layerfilter0PL_newfea_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
:
é
>layerfilter0PL_newfea_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference(layerfilter0PL_newfea_conv_head_0/Conv2D9layerfilter0PL_newfea_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:

Glayerfilter0PL_newfea_conv_head_0/bn/moments/variance/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0

5layerfilter0PL_newfea_conv_head_0/bn/moments/varianceMean>layerfilter0PL_newfea_conv_head_0/bn/moments/SquaredDifferenceGlayerfilter0PL_newfea_conv_head_0/bn/moments/variance/reduction_indices*
T0*&
_output_shapes
:*
	keep_dims(*

Tidx0
°
4layerfilter0PL_newfea_conv_head_0/bn/moments/SqueezeSqueeze1layerfilter0PL_newfea_conv_head_0/bn/moments/mean*
squeeze_dims
 *
T0*
_output_shapes
:
¶
6layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1Squeeze5layerfilter0PL_newfea_conv_head_0/bn/moments/variance*
_output_shapes
:*
squeeze_dims
 *
T0
{
0layerfilter0PL_newfea_conv_head_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

2layerfilter0PL_newfea_conv_head_0/bn/cond/switch_tIdentity2layerfilter0PL_newfea_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter0PL_newfea_conv_head_0/bn/cond/switch_fIdentity0layerfilter0PL_newfea_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter0PL_newfea_conv_head_0/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
Ú
layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes
:
æ
rlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
_output_shapes
:*
shared_name *
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
	container 
Ô
ylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
¤
wlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:*
T0*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
Þ
layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:
ê
tlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:*
shared_name *
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:
Ü
{layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
_output_shapes
:*
use_locking(*
T0*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
ª
ylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:*
T0
Â
Hlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter0PL_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ò
Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter0PL_newfea_conv_head_0/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
¢
Vlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
Ì
Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes
:*
T0
È
_layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id* 
_output_shapes
::*
T0*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
È
alayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze1layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
´
Vlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes
:
Ô
Rlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
_output_shapes
:*
use_locking( *
T0*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
À
Ylayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::*
T0
Ô
Zlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter0PL_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¦
Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
Ò
Zlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes
:*
T0
Î
alayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Î
clayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_11layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
º
Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_output_shapes
:
Ü
Tlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
use_locking( *
T0*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
[layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
ö
Blayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverageNoOpS^layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
©
<layerfilter0PL_newfea_conv_head_0/bn/cond/control_dependencyIdentity2layerfilter0PL_newfea_conv_head_0/bn/cond/switch_tC^layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter0PL_newfea_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter0PL_newfea_conv_head_0/bn/cond/NoOpNoOp3^layerfilter0PL_newfea_conv_head_0/bn/cond/switch_f

>layerfilter0PL_newfea_conv_head_0/bn/cond/control_dependency_1Identity2layerfilter0PL_newfea_conv_head_0/bn/cond/switch_f/^layerfilter0PL_newfea_conv_head_0/bn/cond/NoOp*E
_class;
97loc:@layerfilter0PL_newfea_conv_head_0/bn/cond/switch_f*
_output_shapes
: *
T0

â
/layerfilter0PL_newfea_conv_head_0/bn/cond/MergeMerge>layerfilter0PL_newfea_conv_head_0/bn/cond/control_dependency_1<layerfilter0PL_newfea_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter0PL_newfea_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


4layerfilter0PL_newfea_conv_head_0/bn/cond_1/switch_tIdentity4layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch:1*
_output_shapes
: *
T0


4layerfilter0PL_newfea_conv_head_0/bn/cond_1/switch_fIdentity2layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
o
3layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
Ö
4layerfilter0PL_newfea_conv_head_0/bn/cond_1/IdentityIdentity=layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity/Switch:10^layerfilter0PL_newfea_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
¤
;layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity/SwitchSwitch4layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze3layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
Ú
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity_1Identity?layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:10^layerfilter0PL_newfea_conv_head_0/bn/cond/Merge*
_output_shapes
:*
T0
ª
=layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch6layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_13layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id*I
_class?
=;loc:@layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::*
T0

4layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_1Switchwlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id*
T0*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
£
4layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_2Switchylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::*
T0
Ö
1layerfilter0PL_newfea_conv_head_0/bn/cond_1/MergeMerge4layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_14layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
Ú
3layerfilter0PL_newfea_conv_head_0/bn/cond_1/Merge_1Merge4layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_26layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity_1*
N*
_output_shapes

:: *
T0
y
4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/add/yConst*
_output_shapes
: *
valueB
 *o:*
dtype0
É
2layerfilter0PL_newfea_conv_head_0/bn/batchnorm/addAdd3layerfilter0PL_newfea_conv_head_0/bn/cond_1/Merge_14layerfilter0PL_newfea_conv_head_0/bn/batchnorm/add/y*
_output_shapes
:*
T0

4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/RsqrtRsqrt2layerfilter0PL_newfea_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
:
Å
2layerfilter0PL_newfea_conv_head_0/bn/batchnorm/mulMul4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/Rsqrt/layerfilter0PL_newfea_conv_head_0/bn/gamma/read*
_output_shapes
:*
T0
Ê
4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/mul_1Mul(layerfilter0PL_newfea_conv_head_0/Conv2D2layerfilter0PL_newfea_conv_head_0/bn/batchnorm/mul*&
_output_shapes
:*
T0
Ç
4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/mul_2Mul1layerfilter0PL_newfea_conv_head_0/bn/cond_1/Merge2layerfilter0PL_newfea_conv_head_0/bn/batchnorm/mul*
_output_shapes
:*
T0
Ä
2layerfilter0PL_newfea_conv_head_0/bn/batchnorm/subSub.layerfilter0PL_newfea_conv_head_0/bn/beta/read4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/mul_2*
_output_shapes
:*
T0
Ö
4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/add_1Add4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/mul_12layerfilter0PL_newfea_conv_head_0/bn/batchnorm/sub*&
_output_shapes
:*
T0

&layerfilter0PL_newfea_conv_head_0/ReluRelu4layerfilter0PL_newfea_conv_head_0/bn/batchnorm/add_1*&
_output_shapes
:*
T0
_
	Squeeze_2SqueezeExpandDims_4*
_output_shapes

:*
squeeze_dims
 *
T0
R
ExpandDims_5/dimConst*
value	B : *
dtype0*
_output_shapes
: 
p
ExpandDims_5
ExpandDims	Squeeze_2ExpandDims_5/dim*
T0*"
_output_shapes
:*

Tdim0
M
range/startConst*
value	B : *
dtype0*
_output_shapes
: 
M
range/limitConst*
_output_shapes
: *
value	B :*
dtype0
M
range/deltaConst*
_output_shapes
: *
value	B :*
dtype0
]
rangeRangerange/startrange/limitrange/delta*
_output_shapes
:*

Tidx0
I
mul_3/yConst*
value	B :*
dtype0*
_output_shapes
: 
A
mul_3Mulrangemul_3/y*
T0*
_output_shapes
:
b
Reshape/shapeConst*!
valueB"         *
dtype0*
_output_shapes
:
c
ReshapeReshapemul_3Reshape/shape*
T0*
Tshape0*"
_output_shapes
:
`
Reshape_1/shapeConst*
valueB"ÿÿÿÿ   *
dtype0*
_output_shapes
:
j
	Reshape_1ReshapeExpandDims_5Reshape_1/shape*
T0*
Tshape0*
_output_shapes

:
L
add_5AddTopKV2:1Reshape*
T0*"
_output_shapes
:

O
GatherV2/axisConst*
value	B : *
dtype0*
_output_shapes
: 

GatherV2GatherV2	Reshape_1add_5GatherV2/axis*
Tindices0*
Tparams0*&
_output_shapes
:
*
Taxis0*

batch_dims 
i
Tile_1/multiplesConst*%
valueB"      
      *
dtype0*
_output_shapes
:
q
Tile_1TileExpandDims_4Tile_1/multiples*

Tmultiples0*
T0*&
_output_shapes
:

O
sub_3SubTile_1GatherV2*
T0*&
_output_shapes
:

Ï
Alayerfilter0PL_edgefea_0/weights/Initializer/random_uniform/shapeConst*%
valueB"            *3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
dtype0*
_output_shapes
:
¹
?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/minConst*
valueB
 *Bãè¾*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
dtype0*
_output_shapes
: 
¹
?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/maxConst*
valueB
 *Bãè>*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
dtype0*
_output_shapes
: 
©
Ilayerfilter0PL_edgefea_0/weights/Initializer/random_uniform/RandomUniformRandomUniformAlayerfilter0PL_edgefea_0/weights/Initializer/random_uniform/shape*
T0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
seed2 *
dtype0*&
_output_shapes
:*

seed 

?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/subSub?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/max?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/min*
_output_shapes
: *
T0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights
¸
?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/mulMulIlayerfilter0PL_edgefea_0/weights/Initializer/random_uniform/RandomUniform?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/sub*&
_output_shapes
:*
T0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights
ª
;layerfilter0PL_edgefea_0/weights/Initializer/random_uniformAdd?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/mul?layerfilter0PL_edgefea_0/weights/Initializer/random_uniform/min*
T0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*&
_output_shapes
:
è
 layerfilter0PL_edgefea_0/weights
VariableV2"/device:CPU:0*
shared_name *3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
	container *
shape:*
dtype0*&
_output_shapes
:
®
'layerfilter0PL_edgefea_0/weights/AssignAssign layerfilter0PL_edgefea_0/weights;layerfilter0PL_edgefea_0/weights/Initializer/random_uniform"/device:CPU:0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
validate_shape(*&
_output_shapes
:*
use_locking(*
T0
È
%layerfilter0PL_edgefea_0/weights/readIdentity layerfilter0PL_edgefea_0/weights"/device:CPU:0*
T0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*&
_output_shapes
:
q
layerfilter0PL_edgefea_0/L2LossL2Loss%layerfilter0PL_edgefea_0/weights/read*
_output_shapes
: *
T0
k
&layerfilter0PL_edgefea_0/weight_loss/yConst*
_output_shapes
: *
valueB
 *    *
dtype0

$layerfilter0PL_edgefea_0/weight_lossMullayerfilter0PL_edgefea_0/L2Loss&layerfilter0PL_edgefea_0/weight_loss/y*
T0*
_output_shapes
: 

layerfilter0PL_edgefea_0/Conv2DConv2Dsub_3%layerfilter0PL_edgefea_0/weights/read*
	dilations
*
T0*
data_formatNHWC*
strides
*
explicit_paddings
 *
use_cudnn_on_gpu(*
paddingVALID*&
_output_shapes
:

²
1layerfilter0PL_edgefea_0/biases/Initializer/ConstConst*
valueB*    *2
_class(
&$loc:@layerfilter0PL_edgefea_0/biases*
dtype0*
_output_shapes
:
Î
layerfilter0PL_edgefea_0/biases
VariableV2"/device:CPU:0*2
_class(
&$loc:@layerfilter0PL_edgefea_0/biases*
	container *
shape:*
dtype0*
_output_shapes
:*
shared_name 

&layerfilter0PL_edgefea_0/biases/AssignAssignlayerfilter0PL_edgefea_0/biases1layerfilter0PL_edgefea_0/biases/Initializer/Const"/device:CPU:0*2
_class(
&$loc:@layerfilter0PL_edgefea_0/biases*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
¹
$layerfilter0PL_edgefea_0/biases/readIdentitylayerfilter0PL_edgefea_0/biases"/device:CPU:0*
_output_shapes
:*
T0*2
_class(
&$loc:@layerfilter0PL_edgefea_0/biases
º
 layerfilter0PL_edgefea_0/BiasAddBiasAddlayerfilter0PL_edgefea_0/Conv2D$layerfilter0PL_edgefea_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:

n
!layerfilter0PL_edgefea_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

 layerfilter0PL_edgefea_0/bn/beta
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes
:*
	container 
ù
'layerfilter0PL_edgefea_0/bn/beta/AssignAssign layerfilter0PL_edgefea_0/bn/beta!layerfilter0PL_edgefea_0/bn/Const*
use_locking(*
T0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/bn/beta*
validate_shape(*
_output_shapes
:
­
%layerfilter0PL_edgefea_0/bn/beta/readIdentity layerfilter0PL_edgefea_0/bn/beta*3
_class)
'%loc:@layerfilter0PL_edgefea_0/bn/beta*
_output_shapes
:*
T0
p
#layerfilter0PL_edgefea_0/bn/Const_1Const*
_output_shapes
:*
valueB*  ?*
dtype0

!layerfilter0PL_edgefea_0/bn/gamma
VariableV2*
dtype0*
_output_shapes
:*
	container *
shape:*
shared_name 
þ
(layerfilter0PL_edgefea_0/bn/gamma/AssignAssign!layerfilter0PL_edgefea_0/bn/gamma#layerfilter0PL_edgefea_0/bn/Const_1*
use_locking(*
T0*4
_class*
(&loc:@layerfilter0PL_edgefea_0/bn/gamma*
validate_shape(*
_output_shapes
:
°
&layerfilter0PL_edgefea_0/bn/gamma/readIdentity!layerfilter0PL_edgefea_0/bn/gamma*
T0*4
_class*
(&loc:@layerfilter0PL_edgefea_0/bn/gamma*
_output_shapes
:

:layerfilter0PL_edgefea_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ü
(layerfilter0PL_edgefea_0/bn/moments/meanMean layerfilter0PL_edgefea_0/BiasAdd:layerfilter0PL_edgefea_0/bn/moments/mean/reduction_indices*&
_output_shapes
:*
	keep_dims(*

Tidx0*
T0

0layerfilter0PL_edgefea_0/bn/moments/StopGradientStopGradient(layerfilter0PL_edgefea_0/bn/moments/mean*
T0*&
_output_shapes
:
Ï
5layerfilter0PL_edgefea_0/bn/moments/SquaredDifferenceSquaredDifference layerfilter0PL_edgefea_0/BiasAdd0layerfilter0PL_edgefea_0/bn/moments/StopGradient*
T0*&
_output_shapes
:


>layerfilter0PL_edgefea_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ù
,layerfilter0PL_edgefea_0/bn/moments/varianceMean5layerfilter0PL_edgefea_0/bn/moments/SquaredDifference>layerfilter0PL_edgefea_0/bn/moments/variance/reduction_indices*&
_output_shapes
:*
	keep_dims(*

Tidx0*
T0

+layerfilter0PL_edgefea_0/bn/moments/SqueezeSqueeze(layerfilter0PL_edgefea_0/bn/moments/mean*
_output_shapes
:*
squeeze_dims
 *
T0
¤
-layerfilter0PL_edgefea_0/bn/moments/Squeeze_1Squeeze,layerfilter0PL_edgefea_0/bn/moments/variance*
T0*
_output_shapes
:*
squeeze_dims
 
r
'layerfilter0PL_edgefea_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

)layerfilter0PL_edgefea_0/bn/cond/switch_tIdentity)layerfilter0PL_edgefea_0/bn/cond/Switch:1*
_output_shapes
: *
T0


)layerfilter0PL_edgefea_0/bn/cond/switch_fIdentity'layerfilter0PL_edgefea_0/bn/cond/Switch*
T0
*
_output_shapes
: 
d
(layerfilter0PL_edgefea_0/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
´
rlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes
:
Á
`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:*
shared_name *s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:

glayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragerlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
í
elayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:*
T0
¸
tlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:
Å
blayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes
:*
shared_name 

ilayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragetlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ó
glayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
°
?layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/decayConst*^layerfilter0PL_edgefea_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
À
Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst*^layerfilter0PL_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

Mlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubOlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x?layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
±
Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubXlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Zlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_output_shapes
:

Vlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchelayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read(layerfilter0PL_edgefea_0/bn/cond/pred_id* 
_output_shapes
::*
T0*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage
¤
Xlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch+layerfilter0PL_edgefea_0/bn/moments/Squeeze(layerfilter0PL_edgefea_0/bn/cond/pred_id*
T0*>
_class4
20loc:@layerfilter0PL_edgefea_0/bn/moments/Squeeze* 
_output_shapes
::

Mlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulOlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Mlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes
:
¦
Ilayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubRlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Mlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
use_locking( *
T0*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:

Playerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage(layerfilter0PL_edgefea_0/bn/cond/pred_id*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::*
T0
Â
Qlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst*^layerfilter0PL_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubQlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x?layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
·
Qlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubZlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1\layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes
:*
T0

Xlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchglayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read(layerfilter0PL_edgefea_0/bn/cond/pred_id*
T0*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
ª
Zlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch-layerfilter0PL_edgefea_0/bn/moments/Squeeze_1(layerfilter0PL_edgefea_0/bn/cond/pred_id*@
_class6
42loc:@layerfilter0PL_edgefea_0/bn/moments/Squeeze_1* 
_output_shapes
::*
T0

Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulQlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes
:*
T0
®
Klayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubTlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
_output_shapes
:*
use_locking( *
T0*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage

Rlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage(layerfilter0PL_edgefea_0/bn/cond/pred_id*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::*
T0
Û
9layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverageNoOpJ^layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgL^layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1

3layerfilter0PL_edgefea_0/bn/cond/control_dependencyIdentity)layerfilter0PL_edgefea_0/bn/cond/switch_t:^layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage*
T0
*<
_class2
0.loc:@layerfilter0PL_edgefea_0/bn/cond/switch_t*
_output_shapes
: 
Y
%layerfilter0PL_edgefea_0/bn/cond/NoOpNoOp*^layerfilter0PL_edgefea_0/bn/cond/switch_f
ó
5layerfilter0PL_edgefea_0/bn/cond/control_dependency_1Identity)layerfilter0PL_edgefea_0/bn/cond/switch_f&^layerfilter0PL_edgefea_0/bn/cond/NoOp*
T0
*<
_class2
0.loc:@layerfilter0PL_edgefea_0/bn/cond/switch_f*
_output_shapes
: 
Ç
&layerfilter0PL_edgefea_0/bn/cond/MergeMerge5layerfilter0PL_edgefea_0/bn/cond/control_dependency_13layerfilter0PL_edgefea_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
t
)layerfilter0PL_edgefea_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

+layerfilter0PL_edgefea_0/bn/cond_1/switch_tIdentity+layerfilter0PL_edgefea_0/bn/cond_1/Switch:1*
_output_shapes
: *
T0


+layerfilter0PL_edgefea_0/bn/cond_1/switch_fIdentity)layerfilter0PL_edgefea_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
f
*layerfilter0PL_edgefea_0/bn/cond_1/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

»
+layerfilter0PL_edgefea_0/bn/cond_1/IdentityIdentity4layerfilter0PL_edgefea_0/bn/cond_1/Identity/Switch:1'^layerfilter0PL_edgefea_0/bn/cond/Merge*
_output_shapes
:*
T0

2layerfilter0PL_edgefea_0/bn/cond_1/Identity/SwitchSwitch+layerfilter0PL_edgefea_0/bn/moments/Squeeze*layerfilter0PL_edgefea_0/bn/cond_1/pred_id*
T0*>
_class4
20loc:@layerfilter0PL_edgefea_0/bn/moments/Squeeze* 
_output_shapes
::
¿
-layerfilter0PL_edgefea_0/bn/cond_1/Identity_1Identity6layerfilter0PL_edgefea_0/bn/cond_1/Identity_1/Switch:1'^layerfilter0PL_edgefea_0/bn/cond/Merge*
T0*
_output_shapes
:

4layerfilter0PL_edgefea_0/bn/cond_1/Identity_1/SwitchSwitch-layerfilter0PL_edgefea_0/bn/moments/Squeeze_1*layerfilter0PL_edgefea_0/bn/cond_1/pred_id*@
_class6
42loc:@layerfilter0PL_edgefea_0/bn/moments/Squeeze_1* 
_output_shapes
::*
T0
è
+layerfilter0PL_edgefea_0/bn/cond_1/Switch_1Switchelayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read*layerfilter0PL_edgefea_0/bn/cond_1/pred_id*
T0*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
ì
+layerfilter0PL_edgefea_0/bn/cond_1/Switch_2Switchglayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read*layerfilter0PL_edgefea_0/bn/cond_1/pred_id*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::*
T0
»
(layerfilter0PL_edgefea_0/bn/cond_1/MergeMerge+layerfilter0PL_edgefea_0/bn/cond_1/Switch_1+layerfilter0PL_edgefea_0/bn/cond_1/Identity*
_output_shapes

:: *
T0*
N
¿
*layerfilter0PL_edgefea_0/bn/cond_1/Merge_1Merge+layerfilter0PL_edgefea_0/bn/cond_1/Switch_2-layerfilter0PL_edgefea_0/bn/cond_1/Identity_1*
_output_shapes

:: *
T0*
N
p
+layerfilter0PL_edgefea_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
®
)layerfilter0PL_edgefea_0/bn/batchnorm/addAdd*layerfilter0PL_edgefea_0/bn/cond_1/Merge_1+layerfilter0PL_edgefea_0/bn/batchnorm/add/y*
_output_shapes
:*
T0

+layerfilter0PL_edgefea_0/bn/batchnorm/RsqrtRsqrt)layerfilter0PL_edgefea_0/bn/batchnorm/add*
T0*
_output_shapes
:
ª
)layerfilter0PL_edgefea_0/bn/batchnorm/mulMul+layerfilter0PL_edgefea_0/bn/batchnorm/Rsqrt&layerfilter0PL_edgefea_0/bn/gamma/read*
T0*
_output_shapes
:
°
+layerfilter0PL_edgefea_0/bn/batchnorm/mul_1Mul layerfilter0PL_edgefea_0/BiasAdd)layerfilter0PL_edgefea_0/bn/batchnorm/mul*
T0*&
_output_shapes
:

¬
+layerfilter0PL_edgefea_0/bn/batchnorm/mul_2Mul(layerfilter0PL_edgefea_0/bn/cond_1/Merge)layerfilter0PL_edgefea_0/bn/batchnorm/mul*
_output_shapes
:*
T0
©
)layerfilter0PL_edgefea_0/bn/batchnorm/subSub%layerfilter0PL_edgefea_0/bn/beta/read+layerfilter0PL_edgefea_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:
»
+layerfilter0PL_edgefea_0/bn/batchnorm/add_1Add+layerfilter0PL_edgefea_0/bn/batchnorm/mul_1)layerfilter0PL_edgefea_0/bn/batchnorm/sub*
T0*&
_output_shapes
:


layerfilter0PL_edgefea_0/ReluRelu+layerfilter0PL_edgefea_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:

å
Llayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"            *>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
:
Ï
Jlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *0¿*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ï
Jlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *0?*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ê
Tlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformLlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/shape*
dtype0*&
_output_shapes
:*

seed *
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
seed2 
Ê
Jlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/subSubJlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/maxJlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
_output_shapes
: 
ä
Jlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/mulMulTlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformJlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/sub*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*&
_output_shapes
:*
T0
Ö
Flayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniformAddJlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/mulJlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform/min*&
_output_shapes
:*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights
þ
+layerfilter0PL_self_att_conv_head_0/weights
VariableV2"/device:CPU:0*
dtype0*&
_output_shapes
:*
shared_name *>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
	container *
shape:
Ú
2layerfilter0PL_self_att_conv_head_0/weights/AssignAssign+layerfilter0PL_self_att_conv_head_0/weightsFlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
validate_shape(*&
_output_shapes
:*
use_locking(*
T0
é
0layerfilter0PL_self_att_conv_head_0/weights/readIdentity+layerfilter0PL_self_att_conv_head_0/weights"/device:CPU:0*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*&
_output_shapes
:

*layerfilter0PL_self_att_conv_head_0/L2LossL2Loss0layerfilter0PL_self_att_conv_head_0/weights/read*
_output_shapes
: *
T0
v
1layerfilter0PL_self_att_conv_head_0/weight_loss/yConst*
_output_shapes
: *
valueB
 *    *
dtype0
¶
/layerfilter0PL_self_att_conv_head_0/weight_lossMul*layerfilter0PL_self_att_conv_head_0/L2Loss1layerfilter0PL_self_att_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 
Ç
*layerfilter0PL_self_att_conv_head_0/Conv2DConv2D&layerfilter0PL_newfea_conv_head_0/Relu0layerfilter0PL_self_att_conv_head_0/weights/read*
paddingVALID*&
_output_shapes
:*
	dilations
*
T0*
data_formatNHWC*
strides
*
explicit_paddings
 *
use_cudnn_on_gpu(
È
<layerfilter0PL_self_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *=
_class3
1/loc:@layerfilter0PL_self_att_conv_head_0/biases*
dtype0*
_output_shapes
:
ä
*layerfilter0PL_self_att_conv_head_0/biases
VariableV2"/device:CPU:0*
	container *
shape:*
dtype0*
_output_shapes
:*
shared_name *=
_class3
1/loc:@layerfilter0PL_self_att_conv_head_0/biases
Á
1layerfilter0PL_self_att_conv_head_0/biases/AssignAssign*layerfilter0PL_self_att_conv_head_0/biases<layerfilter0PL_self_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_self_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:
Ú
/layerfilter0PL_self_att_conv_head_0/biases/readIdentity*layerfilter0PL_self_att_conv_head_0/biases"/device:CPU:0*
T0*=
_class3
1/loc:@layerfilter0PL_self_att_conv_head_0/biases*
_output_shapes
:
Û
+layerfilter0PL_self_att_conv_head_0/BiasAddBiasAdd*layerfilter0PL_self_att_conv_head_0/Conv2D/layerfilter0PL_self_att_conv_head_0/biases/read*&
_output_shapes
:*
T0*
data_formatNHWC
y
,layerfilter0PL_self_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

+layerfilter0PL_self_att_conv_head_0/bn/beta
VariableV2*
dtype0*
_output_shapes
:*
	container *
shape:*
shared_name 
¥
2layerfilter0PL_self_att_conv_head_0/bn/beta/AssignAssign+layerfilter0PL_self_att_conv_head_0/bn/beta,layerfilter0PL_self_att_conv_head_0/bn/Const*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
Î
0layerfilter0PL_self_att_conv_head_0/bn/beta/readIdentity+layerfilter0PL_self_att_conv_head_0/bn/beta*
_output_shapes
:*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/bn/beta
{
.layerfilter0PL_self_att_conv_head_0/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

,layerfilter0PL_self_att_conv_head_0/bn/gamma
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes
:*
	container 
ª
3layerfilter0PL_self_att_conv_head_0/bn/gamma/AssignAssign,layerfilter0PL_self_att_conv_head_0/bn/gamma.layerfilter0PL_self_att_conv_head_0/bn/Const_1*
use_locking(*
T0*?
_class5
31loc:@layerfilter0PL_self_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
Ñ
1layerfilter0PL_self_att_conv_head_0/bn/gamma/readIdentity,layerfilter0PL_self_att_conv_head_0/bn/gamma*
T0*?
_class5
31loc:@layerfilter0PL_self_att_conv_head_0/bn/gamma*
_output_shapes
:

Elayerfilter0PL_self_att_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ý
3layerfilter0PL_self_att_conv_head_0/bn/moments/meanMean+layerfilter0PL_self_att_conv_head_0/BiasAddElayerfilter0PL_self_att_conv_head_0/bn/moments/mean/reduction_indices*
T0*&
_output_shapes
:*
	keep_dims(*

Tidx0
±
;layerfilter0PL_self_att_conv_head_0/bn/moments/StopGradientStopGradient3layerfilter0PL_self_att_conv_head_0/bn/moments/mean*&
_output_shapes
:*
T0
ð
@layerfilter0PL_self_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference+layerfilter0PL_self_att_conv_head_0/BiasAdd;layerfilter0PL_self_att_conv_head_0/bn/moments/StopGradient*&
_output_shapes
:*
T0

Ilayerfilter0PL_self_att_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

7layerfilter0PL_self_att_conv_head_0/bn/moments/varianceMean@layerfilter0PL_self_att_conv_head_0/bn/moments/SquaredDifferenceIlayerfilter0PL_self_att_conv_head_0/bn/moments/variance/reduction_indices*
T0*&
_output_shapes
:*
	keep_dims(*

Tidx0
´
6layerfilter0PL_self_att_conv_head_0/bn/moments/SqueezeSqueeze3layerfilter0PL_self_att_conv_head_0/bn/moments/mean*
_output_shapes
:*
squeeze_dims
 *
T0
º
8layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1Squeeze7layerfilter0PL_self_att_conv_head_0/bn/moments/variance*
T0*
_output_shapes
:*
squeeze_dims
 
}
2layerfilter0PL_self_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


4layerfilter0PL_self_att_conv_head_0/bn/cond/switch_tIdentity4layerfilter0PL_self_att_conv_head_0/bn/cond/Switch:1*
_output_shapes
: *
T0


4layerfilter0PL_self_att_conv_head_0/bn/cond/switch_fIdentity2layerfilter0PL_self_att_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
o
3layerfilter0PL_self_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
â
layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes
:
î
vlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes
:*
shared_name 
ä
}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
°
{layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:*
T0
ç
layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:
ó
xlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shared_name *
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes
:
í
layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
_output_shapes
:*
use_locking(*
T0*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
·
}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:*
T0
Æ
Jlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst5^layerfilter0PL_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ö
Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst5^layerfilter0PL_self_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
¨
Xlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubZlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xJlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
Ò
Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subclayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1elayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes
:*
T0
Ô
alayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch{layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Ð
clayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch6layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze3layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id*I
_class?
=;loc:@layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::*
T0
º
Xlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulZlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Xlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
_output_shapes
:*
T0
Þ
Tlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub]layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Xlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
_output_shapes
:*
use_locking( *
T0*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
Ì
[layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage3layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Ø
\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst5^layerfilter0PL_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¬
Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xJlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
Ø
\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subelayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1glayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes
:*
T0
Û
clayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitch}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
elayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch8layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_13layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id*K
_classA
?=loc:@layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::*
T0
À
Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes
:*
T0
ç
Vlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub_layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
_output_shapes
:*
use_locking( *
T0*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
Ó
]layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage3layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::*
T0
ü
Dlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpU^layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgW^layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
±
>layerfilter0PL_self_att_conv_head_0/bn/cond/control_dependencyIdentity4layerfilter0PL_self_att_conv_head_0/bn/cond/switch_tE^layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage*G
_class=
;9loc:@layerfilter0PL_self_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: *
T0

o
0layerfilter0PL_self_att_conv_head_0/bn/cond/NoOpNoOp5^layerfilter0PL_self_att_conv_head_0/bn/cond/switch_f

@layerfilter0PL_self_att_conv_head_0/bn/cond/control_dependency_1Identity4layerfilter0PL_self_att_conv_head_0/bn/cond/switch_f1^layerfilter0PL_self_att_conv_head_0/bn/cond/NoOp*
T0
*G
_class=
;9loc:@layerfilter0PL_self_att_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
è
1layerfilter0PL_self_att_conv_head_0/bn/cond/MergeMerge@layerfilter0PL_self_att_conv_head_0/bn/cond/control_dependency_1>layerfilter0PL_self_att_conv_head_0/bn/cond/control_dependency*
_output_shapes
: : *
T0
*
N

4layerfilter0PL_self_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


6layerfilter0PL_self_att_conv_head_0/bn/cond_1/switch_tIdentity6layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch:1*
_output_shapes
: *
T0


6layerfilter0PL_self_att_conv_head_0/bn/cond_1/switch_fIdentity4layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
q
5layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
Ü
6layerfilter0PL_self_att_conv_head_0/bn/cond_1/IdentityIdentity?layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity/Switch:12^layerfilter0PL_self_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
¬
=layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch6layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze5layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id*I
_class?
=;loc:@layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::*
T0
à
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity_1IdentityAlayerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:12^layerfilter0PL_self_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
²
?layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch8layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_15layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id*K
_classA
?=loc:@layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::*
T0
«
6layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_1Switch{layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read5layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id* 
_output_shapes
::*
T0*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
°
6layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_2Switch}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read5layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ü
3layerfilter0PL_self_att_conv_head_0/bn/cond_1/MergeMerge6layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_16layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity*
_output_shapes

:: *
T0*
N
à
5layerfilter0PL_self_att_conv_head_0/bn/cond_1/Merge_1Merge6layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_28layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity_1*
_output_shapes

:: *
T0*
N
{
6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
Ï
4layerfilter0PL_self_att_conv_head_0/bn/batchnorm/addAdd5layerfilter0PL_self_att_conv_head_0/bn/cond_1/Merge_16layerfilter0PL_self_att_conv_head_0/bn/batchnorm/add/y*
_output_shapes
:*
T0

6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/RsqrtRsqrt4layerfilter0PL_self_att_conv_head_0/bn/batchnorm/add*
_output_shapes
:*
T0
Ë
4layerfilter0PL_self_att_conv_head_0/bn/batchnorm/mulMul6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/Rsqrt1layerfilter0PL_self_att_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:
Ñ
6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/mul_1Mul+layerfilter0PL_self_att_conv_head_0/BiasAdd4layerfilter0PL_self_att_conv_head_0/bn/batchnorm/mul*&
_output_shapes
:*
T0
Í
6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/mul_2Mul3layerfilter0PL_self_att_conv_head_0/bn/cond_1/Merge4layerfilter0PL_self_att_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:
Ê
4layerfilter0PL_self_att_conv_head_0/bn/batchnorm/subSub0layerfilter0PL_self_att_conv_head_0/bn/beta/read6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/mul_2*
_output_shapes
:*
T0
Ü
6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/add_1Add6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/mul_14layerfilter0PL_self_att_conv_head_0/bn/batchnorm/sub*&
_output_shapes
:*
T0

(layerfilter0PL_self_att_conv_head_0/ReluRelu6layerfilter0PL_self_att_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:
å
Llayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"            *>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
dtype0*
_output_shapes
:
Ï
Jlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *0¿*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ï
Jlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
_output_shapes
: *
valueB
 *0?*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
dtype0
Ê
Tlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformLlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/shape*
dtype0*&
_output_shapes
:*

seed *
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
seed2 
Ê
Jlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/subSubJlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/maxJlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
_output_shapes
: 
ä
Jlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/mulMulTlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformJlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/sub*&
_output_shapes
:*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights
Ö
Flayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniformAddJlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/mulJlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*&
_output_shapes
:
þ
+layerfilter0PL_neib_att_conv_head_0/weights
VariableV2"/device:CPU:0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
	container *
shape:*
dtype0*&
_output_shapes
:*
shared_name 
Ú
2layerfilter0PL_neib_att_conv_head_0/weights/AssignAssign+layerfilter0PL_neib_att_conv_head_0/weightsFlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*&
_output_shapes
:*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
validate_shape(
é
0layerfilter0PL_neib_att_conv_head_0/weights/readIdentity+layerfilter0PL_neib_att_conv_head_0/weights"/device:CPU:0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*&
_output_shapes
:*
T0

*layerfilter0PL_neib_att_conv_head_0/L2LossL2Loss0layerfilter0PL_neib_att_conv_head_0/weights/read*
_output_shapes
: *
T0
v
1layerfilter0PL_neib_att_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
¶
/layerfilter0PL_neib_att_conv_head_0/weight_lossMul*layerfilter0PL_neib_att_conv_head_0/L2Loss1layerfilter0PL_neib_att_conv_head_0/weight_loss/y*
_output_shapes
: *
T0
¾
*layerfilter0PL_neib_att_conv_head_0/Conv2DConv2Dlayerfilter0PL_edgefea_0/Relu0layerfilter0PL_neib_att_conv_head_0/weights/read*
data_formatNHWC*
strides
*
explicit_paddings
 *
use_cudnn_on_gpu(*
paddingVALID*&
_output_shapes
:
*
	dilations
*
T0
È
<layerfilter0PL_neib_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *=
_class3
1/loc:@layerfilter0PL_neib_att_conv_head_0/biases*
dtype0*
_output_shapes
:
ä
*layerfilter0PL_neib_att_conv_head_0/biases
VariableV2"/device:CPU:0*
dtype0*
_output_shapes
:*
shared_name *=
_class3
1/loc:@layerfilter0PL_neib_att_conv_head_0/biases*
	container *
shape:
Á
1layerfilter0PL_neib_att_conv_head_0/biases/AssignAssign*layerfilter0PL_neib_att_conv_head_0/biases<layerfilter0PL_neib_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_neib_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:
Ú
/layerfilter0PL_neib_att_conv_head_0/biases/readIdentity*layerfilter0PL_neib_att_conv_head_0/biases"/device:CPU:0*
T0*=
_class3
1/loc:@layerfilter0PL_neib_att_conv_head_0/biases*
_output_shapes
:
Û
+layerfilter0PL_neib_att_conv_head_0/BiasAddBiasAdd*layerfilter0PL_neib_att_conv_head_0/Conv2D/layerfilter0PL_neib_att_conv_head_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:

y
,layerfilter0PL_neib_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

+layerfilter0PL_neib_att_conv_head_0/bn/beta
VariableV2*
dtype0*
_output_shapes
:*
	container *
shape:*
shared_name 
¥
2layerfilter0PL_neib_att_conv_head_0/bn/beta/AssignAssign+layerfilter0PL_neib_att_conv_head_0/bn/beta,layerfilter0PL_neib_att_conv_head_0/bn/Const*
_output_shapes
:*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/bn/beta*
validate_shape(
Î
0layerfilter0PL_neib_att_conv_head_0/bn/beta/readIdentity+layerfilter0PL_neib_att_conv_head_0/bn/beta*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/bn/beta*
_output_shapes
:*
T0
{
.layerfilter0PL_neib_att_conv_head_0/bn/Const_1Const*
_output_shapes
:*
valueB*  ?*
dtype0

,layerfilter0PL_neib_att_conv_head_0/bn/gamma
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes
:*
	container 
ª
3layerfilter0PL_neib_att_conv_head_0/bn/gamma/AssignAssign,layerfilter0PL_neib_att_conv_head_0/bn/gamma.layerfilter0PL_neib_att_conv_head_0/bn/Const_1*
use_locking(*
T0*?
_class5
31loc:@layerfilter0PL_neib_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
Ñ
1layerfilter0PL_neib_att_conv_head_0/bn/gamma/readIdentity,layerfilter0PL_neib_att_conv_head_0/bn/gamma*
_output_shapes
:*
T0*?
_class5
31loc:@layerfilter0PL_neib_att_conv_head_0/bn/gamma

Elayerfilter0PL_neib_att_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ý
3layerfilter0PL_neib_att_conv_head_0/bn/moments/meanMean+layerfilter0PL_neib_att_conv_head_0/BiasAddElayerfilter0PL_neib_att_conv_head_0/bn/moments/mean/reduction_indices*
T0*&
_output_shapes
:*
	keep_dims(*

Tidx0
±
;layerfilter0PL_neib_att_conv_head_0/bn/moments/StopGradientStopGradient3layerfilter0PL_neib_att_conv_head_0/bn/moments/mean*&
_output_shapes
:*
T0
ð
@layerfilter0PL_neib_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference+layerfilter0PL_neib_att_conv_head_0/BiasAdd;layerfilter0PL_neib_att_conv_head_0/bn/moments/StopGradient*&
_output_shapes
:
*
T0

Ilayerfilter0PL_neib_att_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

7layerfilter0PL_neib_att_conv_head_0/bn/moments/varianceMean@layerfilter0PL_neib_att_conv_head_0/bn/moments/SquaredDifferenceIlayerfilter0PL_neib_att_conv_head_0/bn/moments/variance/reduction_indices*&
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
´
6layerfilter0PL_neib_att_conv_head_0/bn/moments/SqueezeSqueeze3layerfilter0PL_neib_att_conv_head_0/bn/moments/mean*
T0*
_output_shapes
:*
squeeze_dims
 
º
8layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1Squeeze7layerfilter0PL_neib_att_conv_head_0/bn/moments/variance*
squeeze_dims
 *
T0*
_output_shapes
:
}
2layerfilter0PL_neib_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


4layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_tIdentity4layerfilter0PL_neib_att_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

4layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_fIdentity2layerfilter0PL_neib_att_conv_head_0/bn/cond/Switch*
_output_shapes
: *
T0

o
3layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

â
layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes
:*
valueB*    *
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0
î
vlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shared_name *
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes
:
ä
}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
°
{layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ç
layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes
:*
valueB*    *
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0
ó
xlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:*
shared_name *
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:
í
layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
·
}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
Jlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst5^layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ö
Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst5^layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¨
Xlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubZlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xJlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
Ò
Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subclayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1elayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_output_shapes
:
Ô
alayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch{layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::*
T0
Ð
clayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch6layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze3layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
º
Xlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulZlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Xlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes
:
Þ
Tlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub]layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Xlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
use_locking( *
T0*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ì
[layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage3layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Ø
\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst5^layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¬
Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xJlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
Ø
\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subelayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1glayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_output_shapes
:
Û
clayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitch}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
elayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch8layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_13layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*K
_classA
?=loc:@layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
À
Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes
:*
T0
ç
Vlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub_layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
_output_shapes
:*
use_locking( *
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
Ó
]layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage3layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
ü
Dlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpU^layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgW^layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
±
>layerfilter0PL_neib_att_conv_head_0/bn/cond/control_dependencyIdentity4layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_tE^layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage*G
_class=
;9loc:@layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: *
T0

o
0layerfilter0PL_neib_att_conv_head_0/bn/cond/NoOpNoOp5^layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_f

@layerfilter0PL_neib_att_conv_head_0/bn/cond/control_dependency_1Identity4layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_f1^layerfilter0PL_neib_att_conv_head_0/bn/cond/NoOp*G
_class=
;9loc:@layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_f*
_output_shapes
: *
T0

è
1layerfilter0PL_neib_att_conv_head_0/bn/cond/MergeMerge@layerfilter0PL_neib_att_conv_head_0/bn/cond/control_dependency_1>layerfilter0PL_neib_att_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 

4layerfilter0PL_neib_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/switch_tIdentity6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch:1*
_output_shapes
: *
T0


6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/switch_fIdentity4layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
q
5layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

Ü
6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/IdentityIdentity?layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity/Switch:12^layerfilter0PL_neib_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
¬
=layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch6layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze5layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
à
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity_1IdentityAlayerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:12^layerfilter0PL_neib_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
²
?layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch8layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_15layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id* 
_output_shapes
::*
T0*K
_classA
?=loc:@layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1
«
6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_1Switch{layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read5layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
°
6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_2Switch}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read5layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id* 
_output_shapes
::*
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
Ü
3layerfilter0PL_neib_att_conv_head_0/bn/cond_1/MergeMerge6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_16layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
à
5layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Merge_1Merge6layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_28layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
{
6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
Ï
4layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/addAdd5layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Merge_16layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/add/y*
_output_shapes
:*
T0

6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/RsqrtRsqrt4layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/add*
_output_shapes
:*
T0
Ë
4layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/mulMul6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/Rsqrt1layerfilter0PL_neib_att_conv_head_0/bn/gamma/read*
_output_shapes
:*
T0
Ñ
6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/mul_1Mul+layerfilter0PL_neib_att_conv_head_0/BiasAdd4layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/mul*&
_output_shapes
:
*
T0
Í
6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/mul_2Mul3layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Merge4layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:
Ê
4layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/subSub0layerfilter0PL_neib_att_conv_head_0/bn/beta/read6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/mul_2*
_output_shapes
:*
T0
Ü
6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/add_1Add6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/mul_14layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:


(layerfilter0PL_neib_att_conv_head_0/ReluRelu6layerfilter0PL_neib_att_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:


add_6Add(layerfilter0PL_self_att_conv_head_0/Relu(layerfilter0PL_neib_att_conv_head_0/Relu*
T0*&
_output_shapes
:

i
transpose_4/permConst*%
valueB"             *
dtype0*
_output_shapes
:
o
transpose_4	Transposeadd_6transpose_4/perm*
T0*&
_output_shapes
:
*
Tperm0
d
	LeakyRelu	LeakyRelutranspose_4*
alpha%ÍÌL>*&
_output_shapes
:
*
T0
N
SoftmaxSoftmax	LeakyRelu*
T0*&
_output_shapes
:


MatMul_1BatchMatMulV2Softmaxlayerfilter0PL_edgefea_0/Relu*&
_output_shapes
:*
adj_x( *
adj_y( *
T0
®
/layerfilter0PL/BiasAdd/biases/Initializer/zerosConst*
valueB*    *0
_class&
$"loc:@layerfilter0PL/BiasAdd/biases*
dtype0*
_output_shapes
:
»
layerfilter0PL/BiasAdd/biases
VariableV2*
dtype0*
_output_shapes
:*
shared_name *0
_class&
$"loc:@layerfilter0PL/BiasAdd/biases*
	container *
shape:
þ
$layerfilter0PL/BiasAdd/biases/AssignAssignlayerfilter0PL/BiasAdd/biases/layerfilter0PL/BiasAdd/biases/Initializer/zeros*
use_locking(*
T0*0
_class&
$"loc:@layerfilter0PL/BiasAdd/biases*
validate_shape(*
_output_shapes
:
¤
"layerfilter0PL/BiasAdd/biases/readIdentitylayerfilter0PL/BiasAdd/biases*0
_class&
$"loc:@layerfilter0PL/BiasAdd/biases*
_output_shapes
:*
T0

layerfilter0PL/BiasAdd/BiasAddBiasAddMatMul_1"layerfilter0PL/BiasAdd/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:
]
ReluRelulayerfilter0PL/BiasAdd/BiasAdd*&
_output_shapes
:*
T0
\
concat/concat_dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
I
concatIdentityRelu*
T0*&
_output_shapes
:
[
ExpandDims_6/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
v
ExpandDims_6
ExpandDimsPlaceholderExpandDims_6/dim*&
_output_shapes
:*

Tdim0*
T0
X
concat_1/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 

concat_1ConcatV2ExpandDims_6concatconcat_1/axis*&
_output_shapes
:*

Tidx0*
T0*
N
^
concat_2/concat_dimConst*
_output_shapes
: *
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0
d
concat_2Identitylayerfilter0PL_edgefea_0/Relu*
T0*&
_output_shapes
:

`
Max/reduction_indicesConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
y
MaxMaxconcat_2Max/reduction_indices*&
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
¯
1gapnet00/weights/Initializer/random_uniform/shapeConst*%
valueB"         @   *#
_class
loc:@gapnet00/weights*
dtype0*
_output_shapes
:

/gapnet00/weights/Initializer/random_uniform/minConst*
valueB
 *R¾*#
_class
loc:@gapnet00/weights*
dtype0*
_output_shapes
: 

/gapnet00/weights/Initializer/random_uniform/maxConst*
valueB
 *R>*#
_class
loc:@gapnet00/weights*
dtype0*
_output_shapes
: 
ù
9gapnet00/weights/Initializer/random_uniform/RandomUniformRandomUniform1gapnet00/weights/Initializer/random_uniform/shape*
seed2 *
dtype0*&
_output_shapes
:@*

seed *
T0*#
_class
loc:@gapnet00/weights
Þ
/gapnet00/weights/Initializer/random_uniform/subSub/gapnet00/weights/Initializer/random_uniform/max/gapnet00/weights/Initializer/random_uniform/min*
_output_shapes
: *
T0*#
_class
loc:@gapnet00/weights
ø
/gapnet00/weights/Initializer/random_uniform/mulMul9gapnet00/weights/Initializer/random_uniform/RandomUniform/gapnet00/weights/Initializer/random_uniform/sub*&
_output_shapes
:@*
T0*#
_class
loc:@gapnet00/weights
ê
+gapnet00/weights/Initializer/random_uniformAdd/gapnet00/weights/Initializer/random_uniform/mul/gapnet00/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet00/weights*&
_output_shapes
:@
È
gapnet00/weights
VariableV2"/device:CPU:0*
dtype0*&
_output_shapes
:@*
shared_name *#
_class
loc:@gapnet00/weights*
	container *
shape:@
î
gapnet00/weights/AssignAssigngapnet00/weights+gapnet00/weights/Initializer/random_uniform"/device:CPU:0*#
_class
loc:@gapnet00/weights*
validate_shape(*&
_output_shapes
:@*
use_locking(*
T0

gapnet00/weights/readIdentitygapnet00/weights"/device:CPU:0*&
_output_shapes
:@*
T0*#
_class
loc:@gapnet00/weights
Q
gapnet00/L2LossL2Lossgapnet00/weights/read*
_output_shapes
: *
T0
[
gapnet00/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
e
gapnet00/weight_lossMulgapnet00/L2Lossgapnet00/weight_loss/y*
_output_shapes
: *
T0
ó
gapnet00/Conv2DConv2Dconcat_1gapnet00/weights/read*
data_formatNHWC*
strides
*
use_cudnn_on_gpu(*
explicit_paddings
 *
paddingVALID*&
_output_shapes
:@*
	dilations
*
T0

!gapnet00/biases/Initializer/ConstConst*
valueB@*    *"
_class
loc:@gapnet00/biases*
dtype0*
_output_shapes
:@
®
gapnet00/biases
VariableV2"/device:CPU:0*
shape:@*
dtype0*
_output_shapes
:@*
shared_name *"
_class
loc:@gapnet00/biases*
	container 
Õ
gapnet00/biases/AssignAssigngapnet00/biases!gapnet00/biases/Initializer/Const"/device:CPU:0*
_output_shapes
:@*
use_locking(*
T0*"
_class
loc:@gapnet00/biases*
validate_shape(

gapnet00/biases/readIdentitygapnet00/biases"/device:CPU:0*
T0*"
_class
loc:@gapnet00/biases*
_output_shapes
:@

gapnet00/BiasAddBiasAddgapnet00/Conv2Dgapnet00/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:@
^
gapnet00/bn/ConstConst*
_output_shapes
:@*
valueB@*    *
dtype0
|
gapnet00/bn/beta
VariableV2*
_output_shapes
:@*
	container *
shape:@*
shared_name *
dtype0
¹
gapnet00/bn/beta/AssignAssigngapnet00/bn/betagapnet00/bn/Const*
_output_shapes
:@*
use_locking(*
T0*#
_class
loc:@gapnet00/bn/beta*
validate_shape(
}
gapnet00/bn/beta/readIdentitygapnet00/bn/beta*
T0*#
_class
loc:@gapnet00/bn/beta*
_output_shapes
:@
`
gapnet00/bn/Const_1Const*
valueB@*  ?*
dtype0*
_output_shapes
:@
}
gapnet00/bn/gamma
VariableV2*
shared_name *
dtype0*
_output_shapes
:@*
	container *
shape:@
¾
gapnet00/bn/gamma/AssignAssigngapnet00/bn/gammagapnet00/bn/Const_1*$
_class
loc:@gapnet00/bn/gamma*
validate_shape(*
_output_shapes
:@*
use_locking(*
T0

gapnet00/bn/gamma/readIdentitygapnet00/bn/gamma*
T0*$
_class
loc:@gapnet00/bn/gamma*
_output_shapes
:@

*gapnet00/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
¬
gapnet00/bn/moments/meanMeangapnet00/BiasAdd*gapnet00/bn/moments/mean/reduction_indices*&
_output_shapes
:@*
	keep_dims(*

Tidx0*
T0
{
 gapnet00/bn/moments/StopGradientStopGradientgapnet00/bn/moments/mean*&
_output_shapes
:@*
T0

%gapnet00/bn/moments/SquaredDifferenceSquaredDifferencegapnet00/BiasAdd gapnet00/bn/moments/StopGradient*&
_output_shapes
:@*
T0

.gapnet00/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
É
gapnet00/bn/moments/varianceMean%gapnet00/bn/moments/SquaredDifference.gapnet00/bn/moments/variance/reduction_indices*&
_output_shapes
:@*
	keep_dims(*

Tidx0*
T0
~
gapnet00/bn/moments/SqueezeSqueezegapnet00/bn/moments/mean*
_output_shapes
:@*
squeeze_dims
 *
T0

gapnet00/bn/moments/Squeeze_1Squeezegapnet00/bn/moments/variance*
_output_shapes
:@*
squeeze_dims
 *
T0
b
gapnet00/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

a
gapnet00/bn/cond/switch_tIdentitygapnet00/bn/cond/Switch:1*
_output_shapes
: *
T0

_
gapnet00/bn/cond/switch_fIdentitygapnet00/bn/cond/Switch*
T0
*
_output_shapes
: 
T
gapnet00/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
ô
Rgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes
:@*
valueB@*    *S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0

@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shared_name *S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:@*
dtype0*
_output_shapes
:@

Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageRgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@

Egapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@*
T0
ø
Tgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB@*    *U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:@

Bgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:@*
shared_name *U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:@

Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageTgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@

Ggapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@*
T0

/gapnet00/bn/cond/ExponentialMovingAverage/decayConst^gapnet00/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
 
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet00/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
×
=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x/gapnet00/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0

?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubHgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes
:@*
T0
±
Fgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchEgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet00/bn/cond/pred_id*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@*
T0
ä
Hgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet00/bn/moments/Squeezegapnet00/bn/cond/pred_id* 
_output_shapes
:@:@*
T0*.
_class$
" loc:@gapnet00/bn/moments/Squeeze
é
=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
_output_shapes
:@*
T0
Ö
9gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubBgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
use_locking( *
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
©
@gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAveragegapnet00/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
¢
Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet00/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Û
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubAgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x/gapnet00/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0

Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubJgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Lgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_output_shapes
:@
·
Hgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchGgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet00/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
ê
Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchgapnet00/bn/moments/Squeeze_1gapnet00/bn/cond/pred_id*
T0*0
_class&
$"loc:@gapnet00/bn/moments/Squeeze_1* 
_output_shapes
:@:@
ï
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulAgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes
:@*
T0
Þ
;gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubDgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@*
use_locking( *
T0
¯
Bgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet00/bn/cond/pred_id*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@*
T0
«
)gapnet00/bn/cond/ExponentialMovingAverageNoOp:^gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg<^gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
Å
#gapnet00/bn/cond/control_dependencyIdentitygapnet00/bn/cond/switch_t*^gapnet00/bn/cond/ExponentialMovingAverage*
T0
*,
_class"
 loc:@gapnet00/bn/cond/switch_t*
_output_shapes
: 
9
gapnet00/bn/cond/NoOpNoOp^gapnet00/bn/cond/switch_f
³
%gapnet00/bn/cond/control_dependency_1Identitygapnet00/bn/cond/switch_f^gapnet00/bn/cond/NoOp*
_output_shapes
: *
T0
*,
_class"
 loc:@gapnet00/bn/cond/switch_f

gapnet00/bn/cond/MergeMerge%gapnet00/bn/cond/control_dependency_1#gapnet00/bn/cond/control_dependency*
_output_shapes
: : *
T0
*
N
d
gapnet00/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

e
gapnet00/bn/cond_1/switch_tIdentitygapnet00/bn/cond_1/Switch:1*
_output_shapes
: *
T0

c
gapnet00/bn/cond_1/switch_fIdentitygapnet00/bn/cond_1/Switch*
_output_shapes
: *
T0

V
gapnet00/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 

gapnet00/bn/cond_1/IdentityIdentity$gapnet00/bn/cond_1/Identity/Switch:1^gapnet00/bn/cond/Merge*
_output_shapes
:@*
T0
À
"gapnet00/bn/cond_1/Identity/SwitchSwitchgapnet00/bn/moments/Squeezegapnet00/bn/cond_1/pred_id*
T0*.
_class$
" loc:@gapnet00/bn/moments/Squeeze* 
_output_shapes
:@:@

gapnet00/bn/cond_1/Identity_1Identity&gapnet00/bn/cond_1/Identity_1/Switch:1^gapnet00/bn/cond/Merge*
_output_shapes
:@*
T0
Æ
$gapnet00/bn/cond_1/Identity_1/SwitchSwitchgapnet00/bn/moments/Squeeze_1gapnet00/bn/cond_1/pred_id*
T0*0
_class&
$"loc:@gapnet00/bn/moments/Squeeze_1* 
_output_shapes
:@:@

gapnet00/bn/cond_1/Switch_1SwitchEgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet00/bn/cond_1/pred_id*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@

gapnet00/bn/cond_1/Switch_2SwitchGgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet00/bn/cond_1/pred_id*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@

gapnet00/bn/cond_1/MergeMergegapnet00/bn/cond_1/Switch_1gapnet00/bn/cond_1/Identity*
T0*
N*
_output_shapes

:@: 

gapnet00/bn/cond_1/Merge_1Mergegapnet00/bn/cond_1/Switch_2gapnet00/bn/cond_1/Identity_1*
N*
_output_shapes

:@: *
T0
`
gapnet00/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
~
gapnet00/bn/batchnorm/addAddgapnet00/bn/cond_1/Merge_1gapnet00/bn/batchnorm/add/y*
T0*
_output_shapes
:@
d
gapnet00/bn/batchnorm/RsqrtRsqrtgapnet00/bn/batchnorm/add*
_output_shapes
:@*
T0
z
gapnet00/bn/batchnorm/mulMulgapnet00/bn/batchnorm/Rsqrtgapnet00/bn/gamma/read*
_output_shapes
:@*
T0

gapnet00/bn/batchnorm/mul_1Mulgapnet00/BiasAddgapnet00/bn/batchnorm/mul*&
_output_shapes
:@*
T0
|
gapnet00/bn/batchnorm/mul_2Mulgapnet00/bn/cond_1/Mergegapnet00/bn/batchnorm/mul*
T0*
_output_shapes
:@
y
gapnet00/bn/batchnorm/subSubgapnet00/bn/beta/readgapnet00/bn/batchnorm/mul_2*
T0*
_output_shapes
:@

gapnet00/bn/batchnorm/add_1Addgapnet00/bn/batchnorm/mul_1gapnet00/bn/batchnorm/sub*&
_output_shapes
:@*
T0
c
gapnet00/ReluRelugapnet00/bn/batchnorm/add_1*&
_output_shapes
:@*
T0
³
3gapnet01PL/weights/Initializer/random_uniform/shapeConst*
_output_shapes
:*%
valueB"      @      *%
_class
loc:@gapnet01PL/weights*
dtype0

1gapnet01PL/weights/Initializer/random_uniform/minConst*
valueB
 *ó5¾*%
_class
loc:@gapnet01PL/weights*
dtype0*
_output_shapes
: 

1gapnet01PL/weights/Initializer/random_uniform/maxConst*
_output_shapes
: *
valueB
 *ó5>*%
_class
loc:@gapnet01PL/weights*
dtype0

;gapnet01PL/weights/Initializer/random_uniform/RandomUniformRandomUniform3gapnet01PL/weights/Initializer/random_uniform/shape*
dtype0*'
_output_shapes
:@*

seed *
T0*%
_class
loc:@gapnet01PL/weights*
seed2 
æ
1gapnet01PL/weights/Initializer/random_uniform/subSub1gapnet01PL/weights/Initializer/random_uniform/max1gapnet01PL/weights/Initializer/random_uniform/min*
T0*%
_class
loc:@gapnet01PL/weights*
_output_shapes
: 

1gapnet01PL/weights/Initializer/random_uniform/mulMul;gapnet01PL/weights/Initializer/random_uniform/RandomUniform1gapnet01PL/weights/Initializer/random_uniform/sub*
T0*%
_class
loc:@gapnet01PL/weights*'
_output_shapes
:@
ó
-gapnet01PL/weights/Initializer/random_uniformAdd1gapnet01PL/weights/Initializer/random_uniform/mul1gapnet01PL/weights/Initializer/random_uniform/min*%
_class
loc:@gapnet01PL/weights*'
_output_shapes
:@*
T0
Î
gapnet01PL/weights
VariableV2"/device:CPU:0*
	container *
shape:@*
dtype0*'
_output_shapes
:@*
shared_name *%
_class
loc:@gapnet01PL/weights
÷
gapnet01PL/weights/AssignAssigngapnet01PL/weights-gapnet01PL/weights/Initializer/random_uniform"/device:CPU:0*
use_locking(*
T0*%
_class
loc:@gapnet01PL/weights*
validate_shape(*'
_output_shapes
:@

gapnet01PL/weights/readIdentitygapnet01PL/weights"/device:CPU:0*'
_output_shapes
:@*
T0*%
_class
loc:@gapnet01PL/weights
U
gapnet01PL/L2LossL2Lossgapnet01PL/weights/read*
T0*
_output_shapes
: 
]
gapnet01PL/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
k
gapnet01PL/weight_lossMulgapnet01PL/L2Lossgapnet01PL/weight_loss/y*
_output_shapes
: *
T0
ý
gapnet01PL/Conv2DConv2Dgapnet00/Relugapnet01PL/weights/read*
paddingVALID*'
_output_shapes
:*
	dilations
*
T0*
data_formatNHWC*
strides
*
explicit_paddings
 *
use_cudnn_on_gpu(

#gapnet01PL/biases/Initializer/ConstConst*
valueB*    *$
_class
loc:@gapnet01PL/biases*
dtype0*
_output_shapes	
:
´
gapnet01PL/biases
VariableV2"/device:CPU:0*
dtype0*
_output_shapes	
:*
shared_name *$
_class
loc:@gapnet01PL/biases*
	container *
shape:
Þ
gapnet01PL/biases/AssignAssigngapnet01PL/biases#gapnet01PL/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*$
_class
loc:@gapnet01PL/biases*
validate_shape(*
_output_shapes	
:

gapnet01PL/biases/readIdentitygapnet01PL/biases"/device:CPU:0*$
_class
loc:@gapnet01PL/biases*
_output_shapes	
:*
T0

gapnet01PL/BiasAddBiasAddgapnet01PL/Conv2Dgapnet01PL/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:
b
gapnet01PL/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:

gapnet01PL/bn/beta
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes	
:*
	container 
Â
gapnet01PL/bn/beta/AssignAssigngapnet01PL/bn/betagapnet01PL/bn/Const*%
_class
loc:@gapnet01PL/bn/beta*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

gapnet01PL/bn/beta/readIdentitygapnet01PL/bn/beta*
T0*%
_class
loc:@gapnet01PL/bn/beta*
_output_shapes	
:
d
gapnet01PL/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:

gapnet01PL/bn/gamma
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes	
:*
	container 
Ç
gapnet01PL/bn/gamma/AssignAssigngapnet01PL/bn/gammagapnet01PL/bn/Const_1*&
_class
loc:@gapnet01PL/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

gapnet01PL/bn/gamma/readIdentitygapnet01PL/bn/gamma*
_output_shapes	
:*
T0*&
_class
loc:@gapnet01PL/bn/gamma

,gapnet01PL/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
³
gapnet01PL/bn/moments/meanMeangapnet01PL/BiasAdd,gapnet01PL/bn/moments/mean/reduction_indices*
	keep_dims(*

Tidx0*
T0*'
_output_shapes
:

"gapnet01PL/bn/moments/StopGradientStopGradientgapnet01PL/bn/moments/mean*
T0*'
_output_shapes
:
¦
'gapnet01PL/bn/moments/SquaredDifferenceSquaredDifferencegapnet01PL/BiasAdd"gapnet01PL/bn/moments/StopGradient*'
_output_shapes
:*
T0

0gapnet01PL/bn/moments/variance/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0
Ð
gapnet01PL/bn/moments/varianceMean'gapnet01PL/bn/moments/SquaredDifference0gapnet01PL/bn/moments/variance/reduction_indices*
T0*'
_output_shapes
:*
	keep_dims(*

Tidx0

gapnet01PL/bn/moments/SqueezeSqueezegapnet01PL/bn/moments/mean*
squeeze_dims
 *
T0*
_output_shapes	
:

gapnet01PL/bn/moments/Squeeze_1Squeezegapnet01PL/bn/moments/variance*
squeeze_dims
 *
T0*
_output_shapes	
:
d
gapnet01PL/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

e
gapnet01PL/bn/cond/switch_tIdentitygapnet01PL/bn/cond/Switch:1*
T0
*
_output_shapes
: 
c
gapnet01PL/bn/cond/switch_fIdentitygapnet01PL/bn/cond/Switch*
_output_shapes
: *
T0

V
gapnet01PL/bn/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

þ
Vgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes	
:

Dgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes	
:*
shared_name *W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:

Kgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverageVgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
_output_shapes	
:*
use_locking(*
T0*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(

Igapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*
T0*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Xgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes	
:*
valueB*    *Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0

Fgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes	
:*
shared_name 
£
Mgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverageXgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
 
Kgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

1gapnet01PL/bn/cond/ExponentialMovingAverage/decayConst^gapnet01PL/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
¤
Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet01PL/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Ý
?gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubAgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x1gapnet01PL/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0

Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubJgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Lgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes	
:*
T0
¿
Hgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchIgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet01PL/bn/cond/pred_id*"
_output_shapes
::*
T0*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage
î
Jgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet01PL/bn/moments/Squeezegapnet01PL/bn/cond/pred_id*
T0*0
_class&
$"loc:@gapnet01PL/bn/moments/Squeeze*"
_output_shapes
::
ð
?gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulAgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1?gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes	
:
á
;gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubDgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1?gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
_output_shapes	
:*
use_locking( *
T0*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage
·
Bgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAveragegapnet01PL/bn/cond/pred_id*
T0*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
¦
Cgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet01PL/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
á
Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubCgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x1gapnet01PL/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 

Cgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubLgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Ngapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes	
:*
T0
Å
Jgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchKgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet01PL/bn/cond/pred_id*
T0*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
ô
Lgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchgapnet01PL/bn/moments/Squeeze_1gapnet01PL/bn/cond/pred_id*"
_output_shapes
::*
T0*2
_class(
&$loc:@gapnet01PL/bn/moments/Squeeze_1
ö
Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulCgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_output_shapes	
:
é
=gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubFgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
use_locking( *
T0*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
½
Dgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet01PL/bn/cond/pred_id*
T0*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
±
+gapnet01PL/bn/cond/ExponentialMovingAverageNoOp<^gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg>^gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
Í
%gapnet01PL/bn/cond/control_dependencyIdentitygapnet01PL/bn/cond/switch_t,^gapnet01PL/bn/cond/ExponentialMovingAverage*
_output_shapes
: *
T0
*.
_class$
" loc:@gapnet01PL/bn/cond/switch_t
=
gapnet01PL/bn/cond/NoOpNoOp^gapnet01PL/bn/cond/switch_f
»
'gapnet01PL/bn/cond/control_dependency_1Identitygapnet01PL/bn/cond/switch_f^gapnet01PL/bn/cond/NoOp*
_output_shapes
: *
T0
*.
_class$
" loc:@gapnet01PL/bn/cond/switch_f

gapnet01PL/bn/cond/MergeMerge'gapnet01PL/bn/cond/control_dependency_1%gapnet01PL/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
f
gapnet01PL/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 
i
gapnet01PL/bn/cond_1/switch_tIdentitygapnet01PL/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
g
gapnet01PL/bn/cond_1/switch_fIdentitygapnet01PL/bn/cond_1/Switch*
T0
*
_output_shapes
: 
X
gapnet01PL/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 

gapnet01PL/bn/cond_1/IdentityIdentity&gapnet01PL/bn/cond_1/Identity/Switch:1^gapnet01PL/bn/cond/Merge*
_output_shapes	
:*
T0
Ê
$gapnet01PL/bn/cond_1/Identity/SwitchSwitchgapnet01PL/bn/moments/Squeezegapnet01PL/bn/cond_1/pred_id*"
_output_shapes
::*
T0*0
_class&
$"loc:@gapnet01PL/bn/moments/Squeeze

gapnet01PL/bn/cond_1/Identity_1Identity(gapnet01PL/bn/cond_1/Identity_1/Switch:1^gapnet01PL/bn/cond/Merge*
T0*
_output_shapes	
:
Ð
&gapnet01PL/bn/cond_1/Identity_1/SwitchSwitchgapnet01PL/bn/moments/Squeeze_1gapnet01PL/bn/cond_1/pred_id*
T0*2
_class(
&$loc:@gapnet01PL/bn/moments/Squeeze_1*"
_output_shapes
::

gapnet01PL/bn/cond_1/Switch_1SwitchIgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet01PL/bn/cond_1/pred_id*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::*
T0

gapnet01PL/bn/cond_1/Switch_2SwitchKgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet01PL/bn/cond_1/pred_id*"
_output_shapes
::*
T0*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage

gapnet01PL/bn/cond_1/MergeMergegapnet01PL/bn/cond_1/Switch_1gapnet01PL/bn/cond_1/Identity*
T0*
N*
_output_shapes
	:: 

gapnet01PL/bn/cond_1/Merge_1Mergegapnet01PL/bn/cond_1/Switch_2gapnet01PL/bn/cond_1/Identity_1*
N*
_output_shapes
	:: *
T0
b
gapnet01PL/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 

gapnet01PL/bn/batchnorm/addAddgapnet01PL/bn/cond_1/Merge_1gapnet01PL/bn/batchnorm/add/y*
T0*
_output_shapes	
:
i
gapnet01PL/bn/batchnorm/RsqrtRsqrtgapnet01PL/bn/batchnorm/add*
_output_shapes	
:*
T0

gapnet01PL/bn/batchnorm/mulMulgapnet01PL/bn/batchnorm/Rsqrtgapnet01PL/bn/gamma/read*
T0*
_output_shapes	
:

gapnet01PL/bn/batchnorm/mul_1Mulgapnet01PL/BiasAddgapnet01PL/bn/batchnorm/mul*
T0*'
_output_shapes
:

gapnet01PL/bn/batchnorm/mul_2Mulgapnet01PL/bn/cond_1/Mergegapnet01PL/bn/batchnorm/mul*
T0*
_output_shapes	
:

gapnet01PL/bn/batchnorm/subSubgapnet01PL/bn/beta/readgapnet01PL/bn/batchnorm/mul_2*
_output_shapes	
:*
T0

gapnet01PL/bn/batchnorm/add_1Addgapnet01PL/bn/batchnorm/mul_1gapnet01PL/bn/batchnorm/sub*
T0*'
_output_shapes
:
h
gapnet01PL/ReluRelugapnet01PL/bn/batchnorm/add_1*'
_output_shapes
:*
T0
b
Max_1/reduction_indicesConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 

Max_1Maxgapnet01PL/ReluMax_1/reduction_indices*'
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
c
	Squeeze_3Squeezegapnet01PL/Relu*
_output_shapes
:	*
squeeze_dims
 *
T0
R
ExpandDims_7/dimConst*
value	B : *
dtype0*
_output_shapes
: 
q
ExpandDims_7
ExpandDims	Squeeze_3ExpandDims_7/dim*
T0*#
_output_shapes
:*

Tdim0
e
transpose_5/permConst*!
valueB"          *
dtype0*
_output_shapes
:
s
transpose_5	TransposeExpandDims_7transpose_5/perm*#
_output_shapes
:*
Tperm0*
T0
{
MatMul_2BatchMatMulV2ExpandDims_7transpose_5*
adj_x( *
adj_y( *
T0*"
_output_shapes
:
L
mul_4/xConst*
valueB
 *   À*
dtype0*
_output_shapes
: 
L
mul_4Mulmul_4/xMatMul_2*"
_output_shapes
:*
T0
N
Square_1SquareExpandDims_7*
T0*#
_output_shapes
:
b
Sum_1/reduction_indicesConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
y
Sum_1SumSquare_1Sum_1/reduction_indices*
	keep_dims(*

Tidx0*
T0*"
_output_shapes
:
e
transpose_6/permConst*!
valueB"          *
dtype0*
_output_shapes
:
k
transpose_6	TransposeSum_1transpose_6/perm*"
_output_shapes
:*
Tperm0*
T0
G
add_7AddSum_1mul_4*"
_output_shapes
:*
T0
M
add_8Addadd_7transpose_6*"
_output_shapes
:*
T0
@
Neg_1Negadd_8*"
_output_shapes
:*
T0
L

TopKV2_1/kConst*
value	B :
*
dtype0*
_output_shapes
: 
n
TopKV2_1TopKV2Neg_1
TopKV2_1/k*
sorted(*
T0*0
_output_shapes
:
:

c
	Squeeze_4Squeezegapnet01PL/Relu*
T0*
_output_shapes
:	*
squeeze_dims
 
R
ExpandDims_8/dimConst*
value	B : *
dtype0*
_output_shapes
: 
q
ExpandDims_8
ExpandDims	Squeeze_4ExpandDims_8/dim*#
_output_shapes
:*

Tdim0*
T0
[
ExpandDims_9/dimConst*
_output_shapes
: *
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0
x
ExpandDims_9
ExpandDimsExpandDims_8ExpandDims_9/dim*

Tdim0*
T0*'
_output_shapes
:
á
Jlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"         @   *<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
dtype0*
_output_shapes
:
Ë
Hlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/minConst*
_output_shapes
: *
valueB
 *ó5¾*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
dtype0
Ë
Hlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *ó5>*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
dtype0*
_output_shapes
: 
Å
Rlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/shape*
dtype0*'
_output_shapes
:@*

seed *
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
seed2 
Â
Hlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/subSubHlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/maxHlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
_output_shapes
: 
Ý
Hlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/mulMulRlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformHlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/sub*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*'
_output_shapes
:@*
T0
Ï
Dlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniformAddHlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/mulHlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*'
_output_shapes
:@
ü
)layerfilter1PL_newfea_conv_head_0/weights
VariableV2"/device:CPU:0*
	container *
shape:@*
dtype0*'
_output_shapes
:@*
shared_name *<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights
Ó
0layerfilter1PL_newfea_conv_head_0/weights/AssignAssign)layerfilter1PL_newfea_conv_head_0/weightsDlayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*'
_output_shapes
:@*
use_locking(*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
validate_shape(
ä
.layerfilter1PL_newfea_conv_head_0/weights/readIdentity)layerfilter1PL_newfea_conv_head_0/weights"/device:CPU:0*'
_output_shapes
:@*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights

(layerfilter1PL_newfea_conv_head_0/L2LossL2Loss.layerfilter1PL_newfea_conv_head_0/weights/read*
_output_shapes
: *
T0
t
/layerfilter1PL_newfea_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter1PL_newfea_conv_head_0/weight_lossMul(layerfilter1PL_newfea_conv_head_0/L2Loss/layerfilter1PL_newfea_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 
©
(layerfilter1PL_newfea_conv_head_0/Conv2DConv2DExpandDims_9.layerfilter1PL_newfea_conv_head_0/weights/read*
data_formatNHWC*
strides
*
use_cudnn_on_gpu(*
explicit_paddings
 *
paddingVALID*&
_output_shapes
:@*
	dilations
*
T0
w
*layerfilter1PL_newfea_conv_head_0/bn/ConstConst*
valueB@*    *
dtype0*
_output_shapes
:@

)layerfilter1PL_newfea_conv_head_0/bn/beta
VariableV2*
shape:@*
shared_name *
dtype0*
_output_shapes
:@*
	container 

0layerfilter1PL_newfea_conv_head_0/bn/beta/AssignAssign)layerfilter1PL_newfea_conv_head_0/bn/beta*layerfilter1PL_newfea_conv_head_0/bn/Const*
use_locking(*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:@
È
.layerfilter1PL_newfea_conv_head_0/bn/beta/readIdentity)layerfilter1PL_newfea_conv_head_0/bn/beta*
_output_shapes
:@*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/bn/beta
y
,layerfilter1PL_newfea_conv_head_0/bn/Const_1Const*
valueB@*  ?*
dtype0*
_output_shapes
:@

*layerfilter1PL_newfea_conv_head_0/bn/gamma
VariableV2*
shared_name *
dtype0*
_output_shapes
:@*
	container *
shape:@
¢
1layerfilter1PL_newfea_conv_head_0/bn/gamma/AssignAssign*layerfilter1PL_newfea_conv_head_0/bn/gamma,layerfilter1PL_newfea_conv_head_0/bn/Const_1*
_output_shapes
:@*
use_locking(*
T0*=
_class3
1/loc:@layerfilter1PL_newfea_conv_head_0/bn/gamma*
validate_shape(
Ë
/layerfilter1PL_newfea_conv_head_0/bn/gamma/readIdentity*layerfilter1PL_newfea_conv_head_0/bn/gamma*
T0*=
_class3
1/loc:@layerfilter1PL_newfea_conv_head_0/bn/gamma*
_output_shapes
:@

Clayerfilter1PL_newfea_conv_head_0/bn/moments/mean/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0
ö
1layerfilter1PL_newfea_conv_head_0/bn/moments/meanMean(layerfilter1PL_newfea_conv_head_0/Conv2DClayerfilter1PL_newfea_conv_head_0/bn/moments/mean/reduction_indices*&
_output_shapes
:@*
	keep_dims(*

Tidx0*
T0
­
9layerfilter1PL_newfea_conv_head_0/bn/moments/StopGradientStopGradient1layerfilter1PL_newfea_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
:@
é
>layerfilter1PL_newfea_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference(layerfilter1PL_newfea_conv_head_0/Conv2D9layerfilter1PL_newfea_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:@

Glayerfilter1PL_newfea_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

5layerfilter1PL_newfea_conv_head_0/bn/moments/varianceMean>layerfilter1PL_newfea_conv_head_0/bn/moments/SquaredDifferenceGlayerfilter1PL_newfea_conv_head_0/bn/moments/variance/reduction_indices*&
_output_shapes
:@*
	keep_dims(*

Tidx0*
T0
°
4layerfilter1PL_newfea_conv_head_0/bn/moments/SqueezeSqueeze1layerfilter1PL_newfea_conv_head_0/bn/moments/mean*
squeeze_dims
 *
T0*
_output_shapes
:@
¶
6layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1Squeeze5layerfilter1PL_newfea_conv_head_0/bn/moments/variance*
_output_shapes
:@*
squeeze_dims
 *
T0
{
0layerfilter1PL_newfea_conv_head_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


2layerfilter1PL_newfea_conv_head_0/bn/cond/switch_tIdentity2layerfilter1PL_newfea_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter1PL_newfea_conv_head_0/bn/cond/switch_fIdentity0layerfilter1PL_newfea_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter1PL_newfea_conv_head_0/bn/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

Ú
layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes
:@*
valueB@*    *
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0
æ
rlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:@*
shared_name *
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:@
Ô
ylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
¤
wlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@*
T0
Þ
layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB@*    *
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:@
ê
tlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:@*
shared_name *
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:@
Ü
{layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
_output_shapes
:@*
use_locking(*
T0*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
ª
ylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@*
T0*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
Â
Hlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter1PL_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ò
Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter1PL_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¢
Vlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
Ì
Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_output_shapes
:@
È
_layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
È
alayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze1layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id*G
_class=
;9loc:@layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
:@:@*
T0
´
Vlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
_output_shapes
:@*
T0
Ô
Rlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
use_locking( *
T0*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
À
Ylayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
Ô
Zlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter1PL_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¦
Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
Ò
Zlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes
:@*
T0
Î
alayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Î
clayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_11layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id*I
_class?
=;loc:@layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
:@:@*
T0
º
Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes
:@*
T0
Ü
Tlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@*
use_locking( *
T0
Æ
[layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id* 
_output_shapes
:@:@*
T0*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
ö
Blayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverageNoOpS^layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
©
<layerfilter1PL_newfea_conv_head_0/bn/cond/control_dependencyIdentity2layerfilter1PL_newfea_conv_head_0/bn/cond/switch_tC^layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter1PL_newfea_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter1PL_newfea_conv_head_0/bn/cond/NoOpNoOp3^layerfilter1PL_newfea_conv_head_0/bn/cond/switch_f

>layerfilter1PL_newfea_conv_head_0/bn/cond/control_dependency_1Identity2layerfilter1PL_newfea_conv_head_0/bn/cond/switch_f/^layerfilter1PL_newfea_conv_head_0/bn/cond/NoOp*
T0
*E
_class;
97loc:@layerfilter1PL_newfea_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
â
/layerfilter1PL_newfea_conv_head_0/bn/cond/MergeMerge>layerfilter1PL_newfea_conv_head_0/bn/cond/control_dependency_1<layerfilter1PL_newfea_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter1PL_newfea_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


4layerfilter1PL_newfea_conv_head_0/bn/cond_1/switch_tIdentity4layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch:1*
_output_shapes
: *
T0


4layerfilter1PL_newfea_conv_head_0/bn/cond_1/switch_fIdentity2layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch*
_output_shapes
: *
T0

o
3layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

Ö
4layerfilter1PL_newfea_conv_head_0/bn/cond_1/IdentityIdentity=layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity/Switch:10^layerfilter1PL_newfea_conv_head_0/bn/cond/Merge*
_output_shapes
:@*
T0
¤
;layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity/SwitchSwitch4layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze3layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
:@:@
Ú
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity_1Identity?layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:10^layerfilter1PL_newfea_conv_head_0/bn/cond/Merge*
_output_shapes
:@*
T0
ª
=layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch6layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_13layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id*I
_class?
=;loc:@layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
:@:@*
T0

4layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_1Switchwlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@*
T0
£
4layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_2Switchylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id*
T0*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Ö
1layerfilter1PL_newfea_conv_head_0/bn/cond_1/MergeMerge4layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_14layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity*
N*
_output_shapes

:@: *
T0
Ú
3layerfilter1PL_newfea_conv_head_0/bn/cond_1/Merge_1Merge4layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_26layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:@: 
y
4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
É
2layerfilter1PL_newfea_conv_head_0/bn/batchnorm/addAdd3layerfilter1PL_newfea_conv_head_0/bn/cond_1/Merge_14layerfilter1PL_newfea_conv_head_0/bn/batchnorm/add/y*
_output_shapes
:@*
T0

4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/RsqrtRsqrt2layerfilter1PL_newfea_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
:@
Å
2layerfilter1PL_newfea_conv_head_0/bn/batchnorm/mulMul4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/Rsqrt/layerfilter1PL_newfea_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:@
Ê
4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/mul_1Mul(layerfilter1PL_newfea_conv_head_0/Conv2D2layerfilter1PL_newfea_conv_head_0/bn/batchnorm/mul*&
_output_shapes
:@*
T0
Ç
4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/mul_2Mul1layerfilter1PL_newfea_conv_head_0/bn/cond_1/Merge2layerfilter1PL_newfea_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:@
Ä
2layerfilter1PL_newfea_conv_head_0/bn/batchnorm/subSub.layerfilter1PL_newfea_conv_head_0/bn/beta/read4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/mul_2*
_output_shapes
:@*
T0
Ö
4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/add_1Add4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/mul_12layerfilter1PL_newfea_conv_head_0/bn/batchnorm/sub*&
_output_shapes
:@*
T0

&layerfilter1PL_newfea_conv_head_0/ReluRelu4layerfilter1PL_newfea_conv_head_0/bn/batchnorm/add_1*&
_output_shapes
:@*
T0
`
	Squeeze_5SqueezeExpandDims_9*
T0*
_output_shapes
:	*
squeeze_dims
 
S
ExpandDims_10/dimConst*
value	B : *
dtype0*
_output_shapes
: 
s
ExpandDims_10
ExpandDims	Squeeze_5ExpandDims_10/dim*

Tdim0*
T0*#
_output_shapes
:
O
range_1/startConst*
value	B : *
dtype0*
_output_shapes
: 
O
range_1/limitConst*
value	B :*
dtype0*
_output_shapes
: 
O
range_1/deltaConst*
value	B :*
dtype0*
_output_shapes
: 
e
range_1Rangerange_1/startrange_1/limitrange_1/delta*

Tidx0*
_output_shapes
:
I
mul_5/yConst*
value	B :*
dtype0*
_output_shapes
: 
C
mul_5Mulrange_1mul_5/y*
T0*
_output_shapes
:
d
Reshape_2/shapeConst*
_output_shapes
:*!
valueB"         *
dtype0
g
	Reshape_2Reshapemul_5Reshape_2/shape*
T0*
Tshape0*"
_output_shapes
:
`
Reshape_3/shapeConst*
_output_shapes
:*
valueB"ÿÿÿÿ   *
dtype0
l
	Reshape_3ReshapeExpandDims_10Reshape_3/shape*
_output_shapes
:	*
T0*
Tshape0
P
add_9Add
TopKV2_1:1	Reshape_2*
T0*"
_output_shapes
:

Q
GatherV2_1/axisConst*
_output_shapes
: *
value	B : *
dtype0
 

GatherV2_1GatherV2	Reshape_3add_9GatherV2_1/axis*
Taxis0*

batch_dims *
Tindices0*
Tparams0*'
_output_shapes
:

i
Tile_2/multiplesConst*
_output_shapes
:*%
valueB"      
      *
dtype0
r
Tile_2TileExpandDims_9Tile_2/multiples*'
_output_shapes
:
*

Tmultiples0*
T0
R
sub_4SubTile_2
GatherV2_1*'
_output_shapes
:
*
T0
Ï
Alayerfilter1PL_edgefea_0/weights/Initializer/random_uniform/shapeConst*%
valueB"         @   *3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
dtype0*
_output_shapes
:
¹
?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/minConst*
valueB
 *ó5¾*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
dtype0*
_output_shapes
: 
¹
?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/maxConst*
valueB
 *ó5>*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
dtype0*
_output_shapes
: 
ª
Ilayerfilter1PL_edgefea_0/weights/Initializer/random_uniform/RandomUniformRandomUniformAlayerfilter1PL_edgefea_0/weights/Initializer/random_uniform/shape*'
_output_shapes
:@*

seed *
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
seed2 *
dtype0

?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/subSub?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/max?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/min*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
_output_shapes
: *
T0
¹
?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/mulMulIlayerfilter1PL_edgefea_0/weights/Initializer/random_uniform/RandomUniform?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/sub*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*'
_output_shapes
:@*
T0
«
;layerfilter1PL_edgefea_0/weights/Initializer/random_uniformAdd?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/mul?layerfilter1PL_edgefea_0/weights/Initializer/random_uniform/min*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*'
_output_shapes
:@*
T0
ê
 layerfilter1PL_edgefea_0/weights
VariableV2"/device:CPU:0*
dtype0*'
_output_shapes
:@*
shared_name *3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
	container *
shape:@
¯
'layerfilter1PL_edgefea_0/weights/AssignAssign layerfilter1PL_edgefea_0/weights;layerfilter1PL_edgefea_0/weights/Initializer/random_uniform"/device:CPU:0*
use_locking(*
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
validate_shape(*'
_output_shapes
:@
É
%layerfilter1PL_edgefea_0/weights/readIdentity layerfilter1PL_edgefea_0/weights"/device:CPU:0*'
_output_shapes
:@*
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights
q
layerfilter1PL_edgefea_0/L2LossL2Loss%layerfilter1PL_edgefea_0/weights/read*
T0*
_output_shapes
: 
k
&layerfilter1PL_edgefea_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 

$layerfilter1PL_edgefea_0/weight_lossMullayerfilter1PL_edgefea_0/L2Loss&layerfilter1PL_edgefea_0/weight_loss/y*
_output_shapes
: *
T0

layerfilter1PL_edgefea_0/Conv2DConv2Dsub_4%layerfilter1PL_edgefea_0/weights/read*
use_cudnn_on_gpu(*
explicit_paddings
 *
paddingVALID*&
_output_shapes
:
@*
	dilations
*
T0*
strides
*
data_formatNHWC
²
1layerfilter1PL_edgefea_0/biases/Initializer/ConstConst*
valueB@*    *2
_class(
&$loc:@layerfilter1PL_edgefea_0/biases*
dtype0*
_output_shapes
:@
Î
layerfilter1PL_edgefea_0/biases
VariableV2"/device:CPU:0*
dtype0*
_output_shapes
:@*
shared_name *2
_class(
&$loc:@layerfilter1PL_edgefea_0/biases*
	container *
shape:@

&layerfilter1PL_edgefea_0/biases/AssignAssignlayerfilter1PL_edgefea_0/biases1layerfilter1PL_edgefea_0/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*2
_class(
&$loc:@layerfilter1PL_edgefea_0/biases*
validate_shape(*
_output_shapes
:@
¹
$layerfilter1PL_edgefea_0/biases/readIdentitylayerfilter1PL_edgefea_0/biases"/device:CPU:0*
_output_shapes
:@*
T0*2
_class(
&$loc:@layerfilter1PL_edgefea_0/biases
º
 layerfilter1PL_edgefea_0/BiasAddBiasAddlayerfilter1PL_edgefea_0/Conv2D$layerfilter1PL_edgefea_0/biases/read*&
_output_shapes
:
@*
T0*
data_formatNHWC
n
!layerfilter1PL_edgefea_0/bn/ConstConst*
_output_shapes
:@*
valueB@*    *
dtype0

 layerfilter1PL_edgefea_0/bn/beta
VariableV2*
shape:@*
shared_name *
dtype0*
_output_shapes
:@*
	container 
ù
'layerfilter1PL_edgefea_0/bn/beta/AssignAssign layerfilter1PL_edgefea_0/bn/beta!layerfilter1PL_edgefea_0/bn/Const*
_output_shapes
:@*
use_locking(*
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/bn/beta*
validate_shape(
­
%layerfilter1PL_edgefea_0/bn/beta/readIdentity layerfilter1PL_edgefea_0/bn/beta*3
_class)
'%loc:@layerfilter1PL_edgefea_0/bn/beta*
_output_shapes
:@*
T0
p
#layerfilter1PL_edgefea_0/bn/Const_1Const*
_output_shapes
:@*
valueB@*  ?*
dtype0

!layerfilter1PL_edgefea_0/bn/gamma
VariableV2*
dtype0*
_output_shapes
:@*
	container *
shape:@*
shared_name 
þ
(layerfilter1PL_edgefea_0/bn/gamma/AssignAssign!layerfilter1PL_edgefea_0/bn/gamma#layerfilter1PL_edgefea_0/bn/Const_1*
_output_shapes
:@*
use_locking(*
T0*4
_class*
(&loc:@layerfilter1PL_edgefea_0/bn/gamma*
validate_shape(
°
&layerfilter1PL_edgefea_0/bn/gamma/readIdentity!layerfilter1PL_edgefea_0/bn/gamma*4
_class*
(&loc:@layerfilter1PL_edgefea_0/bn/gamma*
_output_shapes
:@*
T0

:layerfilter1PL_edgefea_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ü
(layerfilter1PL_edgefea_0/bn/moments/meanMean layerfilter1PL_edgefea_0/BiasAdd:layerfilter1PL_edgefea_0/bn/moments/mean/reduction_indices*&
_output_shapes
:@*
	keep_dims(*

Tidx0*
T0

0layerfilter1PL_edgefea_0/bn/moments/StopGradientStopGradient(layerfilter1PL_edgefea_0/bn/moments/mean*
T0*&
_output_shapes
:@
Ï
5layerfilter1PL_edgefea_0/bn/moments/SquaredDifferenceSquaredDifference layerfilter1PL_edgefea_0/BiasAdd0layerfilter1PL_edgefea_0/bn/moments/StopGradient*
T0*&
_output_shapes
:
@

>layerfilter1PL_edgefea_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ù
,layerfilter1PL_edgefea_0/bn/moments/varianceMean5layerfilter1PL_edgefea_0/bn/moments/SquaredDifference>layerfilter1PL_edgefea_0/bn/moments/variance/reduction_indices*
	keep_dims(*

Tidx0*
T0*&
_output_shapes
:@

+layerfilter1PL_edgefea_0/bn/moments/SqueezeSqueeze(layerfilter1PL_edgefea_0/bn/moments/mean*
T0*
_output_shapes
:@*
squeeze_dims
 
¤
-layerfilter1PL_edgefea_0/bn/moments/Squeeze_1Squeeze,layerfilter1PL_edgefea_0/bn/moments/variance*
_output_shapes
:@*
squeeze_dims
 *
T0
r
'layerfilter1PL_edgefea_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


)layerfilter1PL_edgefea_0/bn/cond/switch_tIdentity)layerfilter1PL_edgefea_0/bn/cond/Switch:1*
_output_shapes
: *
T0


)layerfilter1PL_edgefea_0/bn/cond/switch_fIdentity'layerfilter1PL_edgefea_0/bn/cond/Switch*
T0
*
_output_shapes
: 
d
(layerfilter1PL_edgefea_0/bn/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

´
rlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB@*    *s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes
:@
Á
`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
	container *
shape:@*
dtype0*
_output_shapes
:@*
shared_name *s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage

glayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragerlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
í
elayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@*
T0
¸
tlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB@*    *u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:@
Å
blayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shared_name *u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:@*
dtype0*
_output_shapes
:@

ilayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragetlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@*
use_locking(*
T0
ó
glayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
°
?layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/decayConst*^layerfilter1PL_edgefea_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
À
Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst*^layerfilter1PL_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

Mlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubOlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x?layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
±
Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubXlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Zlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_output_shapes
:@

Vlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchelayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read(layerfilter1PL_edgefea_0/bn/cond/pred_id*
T0*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
¤
Xlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch+layerfilter1PL_edgefea_0/bn/moments/Squeeze(layerfilter1PL_edgefea_0/bn/cond/pred_id*
T0*>
_class4
20loc:@layerfilter1PL_edgefea_0/bn/moments/Squeeze* 
_output_shapes
:@:@

Mlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulOlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Mlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
_output_shapes
:@*
T0
¦
Ilayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubRlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Mlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
use_locking( *
T0*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

Playerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage(layerfilter1PL_edgefea_0/bn/cond/pred_id*
T0*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
Â
Qlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst*^layerfilter1PL_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 

Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubQlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x?layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
·
Qlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubZlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1\layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_output_shapes
:@

Xlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchglayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read(layerfilter1PL_edgefea_0/bn/cond/pred_id*
T0*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
ª
Zlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch-layerfilter1PL_edgefea_0/bn/moments/Squeeze_1(layerfilter1PL_edgefea_0/bn/cond/pred_id*
T0*@
_class6
42loc:@layerfilter1PL_edgefea_0/bn/moments/Squeeze_1* 
_output_shapes
:@:@

Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulQlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_output_shapes
:@
®
Klayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubTlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
use_locking( *
T0*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Rlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage(layerfilter1PL_edgefea_0/bn/cond/pred_id*
T0*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Û
9layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverageNoOpJ^layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgL^layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1

3layerfilter1PL_edgefea_0/bn/cond/control_dependencyIdentity)layerfilter1PL_edgefea_0/bn/cond/switch_t:^layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage*
T0
*<
_class2
0.loc:@layerfilter1PL_edgefea_0/bn/cond/switch_t*
_output_shapes
: 
Y
%layerfilter1PL_edgefea_0/bn/cond/NoOpNoOp*^layerfilter1PL_edgefea_0/bn/cond/switch_f
ó
5layerfilter1PL_edgefea_0/bn/cond/control_dependency_1Identity)layerfilter1PL_edgefea_0/bn/cond/switch_f&^layerfilter1PL_edgefea_0/bn/cond/NoOp*
_output_shapes
: *
T0
*<
_class2
0.loc:@layerfilter1PL_edgefea_0/bn/cond/switch_f
Ç
&layerfilter1PL_edgefea_0/bn/cond/MergeMerge5layerfilter1PL_edgefea_0/bn/cond/control_dependency_13layerfilter1PL_edgefea_0/bn/cond/control_dependency*
_output_shapes
: : *
T0
*
N
t
)layerfilter1PL_edgefea_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

+layerfilter1PL_edgefea_0/bn/cond_1/switch_tIdentity+layerfilter1PL_edgefea_0/bn/cond_1/Switch:1*
_output_shapes
: *
T0


+layerfilter1PL_edgefea_0/bn/cond_1/switch_fIdentity)layerfilter1PL_edgefea_0/bn/cond_1/Switch*
_output_shapes
: *
T0

f
*layerfilter1PL_edgefea_0/bn/cond_1/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

»
+layerfilter1PL_edgefea_0/bn/cond_1/IdentityIdentity4layerfilter1PL_edgefea_0/bn/cond_1/Identity/Switch:1'^layerfilter1PL_edgefea_0/bn/cond/Merge*
_output_shapes
:@*
T0

2layerfilter1PL_edgefea_0/bn/cond_1/Identity/SwitchSwitch+layerfilter1PL_edgefea_0/bn/moments/Squeeze*layerfilter1PL_edgefea_0/bn/cond_1/pred_id*>
_class4
20loc:@layerfilter1PL_edgefea_0/bn/moments/Squeeze* 
_output_shapes
:@:@*
T0
¿
-layerfilter1PL_edgefea_0/bn/cond_1/Identity_1Identity6layerfilter1PL_edgefea_0/bn/cond_1/Identity_1/Switch:1'^layerfilter1PL_edgefea_0/bn/cond/Merge*
_output_shapes
:@*
T0

4layerfilter1PL_edgefea_0/bn/cond_1/Identity_1/SwitchSwitch-layerfilter1PL_edgefea_0/bn/moments/Squeeze_1*layerfilter1PL_edgefea_0/bn/cond_1/pred_id* 
_output_shapes
:@:@*
T0*@
_class6
42loc:@layerfilter1PL_edgefea_0/bn/moments/Squeeze_1
è
+layerfilter1PL_edgefea_0/bn/cond_1/Switch_1Switchelayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read*layerfilter1PL_edgefea_0/bn/cond_1/pred_id*
T0*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
ì
+layerfilter1PL_edgefea_0/bn/cond_1/Switch_2Switchglayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read*layerfilter1PL_edgefea_0/bn/cond_1/pred_id*
T0*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
»
(layerfilter1PL_edgefea_0/bn/cond_1/MergeMerge+layerfilter1PL_edgefea_0/bn/cond_1/Switch_1+layerfilter1PL_edgefea_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:@: 
¿
*layerfilter1PL_edgefea_0/bn/cond_1/Merge_1Merge+layerfilter1PL_edgefea_0/bn/cond_1/Switch_2-layerfilter1PL_edgefea_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:@: 
p
+layerfilter1PL_edgefea_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
®
)layerfilter1PL_edgefea_0/bn/batchnorm/addAdd*layerfilter1PL_edgefea_0/bn/cond_1/Merge_1+layerfilter1PL_edgefea_0/bn/batchnorm/add/y*
T0*
_output_shapes
:@

+layerfilter1PL_edgefea_0/bn/batchnorm/RsqrtRsqrt)layerfilter1PL_edgefea_0/bn/batchnorm/add*
_output_shapes
:@*
T0
ª
)layerfilter1PL_edgefea_0/bn/batchnorm/mulMul+layerfilter1PL_edgefea_0/bn/batchnorm/Rsqrt&layerfilter1PL_edgefea_0/bn/gamma/read*
_output_shapes
:@*
T0
°
+layerfilter1PL_edgefea_0/bn/batchnorm/mul_1Mul layerfilter1PL_edgefea_0/BiasAdd)layerfilter1PL_edgefea_0/bn/batchnorm/mul*&
_output_shapes
:
@*
T0
¬
+layerfilter1PL_edgefea_0/bn/batchnorm/mul_2Mul(layerfilter1PL_edgefea_0/bn/cond_1/Merge)layerfilter1PL_edgefea_0/bn/batchnorm/mul*
T0*
_output_shapes
:@
©
)layerfilter1PL_edgefea_0/bn/batchnorm/subSub%layerfilter1PL_edgefea_0/bn/beta/read+layerfilter1PL_edgefea_0/bn/batchnorm/mul_2*
_output_shapes
:@*
T0
»
+layerfilter1PL_edgefea_0/bn/batchnorm/add_1Add+layerfilter1PL_edgefea_0/bn/batchnorm/mul_1)layerfilter1PL_edgefea_0/bn/batchnorm/sub*&
_output_shapes
:
@*
T0

layerfilter1PL_edgefea_0/ReluRelu+layerfilter1PL_edgefea_0/bn/batchnorm/add_1*&
_output_shapes
:
@*
T0
å
Llayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"      @      *>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
:
Ï
Jlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *¾*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ï
Jlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *>*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ê
Tlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformLlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/shape*
dtype0*&
_output_shapes
:@*

seed *
T0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
seed2 
Ê
Jlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/subSubJlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/maxJlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
_output_shapes
: 
ä
Jlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/mulMulTlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformJlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/sub*&
_output_shapes
:@*
T0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights
Ö
Flayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniformAddJlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/mulJlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform/min*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*&
_output_shapes
:@*
T0
þ
+layerfilter1PL_self_att_conv_head_0/weights
VariableV2"/device:CPU:0*
shared_name *>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
	container *
shape:@*
dtype0*&
_output_shapes
:@
Ú
2layerfilter1PL_self_att_conv_head_0/weights/AssignAssign+layerfilter1PL_self_att_conv_head_0/weightsFlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*&
_output_shapes
:@*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
validate_shape(
é
0layerfilter1PL_self_att_conv_head_0/weights/readIdentity+layerfilter1PL_self_att_conv_head_0/weights"/device:CPU:0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*&
_output_shapes
:@*
T0

*layerfilter1PL_self_att_conv_head_0/L2LossL2Loss0layerfilter1PL_self_att_conv_head_0/weights/read*
T0*
_output_shapes
: 
v
1layerfilter1PL_self_att_conv_head_0/weight_loss/yConst*
_output_shapes
: *
valueB
 *    *
dtype0
¶
/layerfilter1PL_self_att_conv_head_0/weight_lossMul*layerfilter1PL_self_att_conv_head_0/L2Loss1layerfilter1PL_self_att_conv_head_0/weight_loss/y*
_output_shapes
: *
T0
Ç
*layerfilter1PL_self_att_conv_head_0/Conv2DConv2D&layerfilter1PL_newfea_conv_head_0/Relu0layerfilter1PL_self_att_conv_head_0/weights/read*&
_output_shapes
:*
	dilations
*
T0*
strides
*
data_formatNHWC*
explicit_paddings
 *
use_cudnn_on_gpu(*
paddingVALID
È
<layerfilter1PL_self_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *=
_class3
1/loc:@layerfilter1PL_self_att_conv_head_0/biases*
dtype0*
_output_shapes
:
ä
*layerfilter1PL_self_att_conv_head_0/biases
VariableV2"/device:CPU:0*
shared_name *=
_class3
1/loc:@layerfilter1PL_self_att_conv_head_0/biases*
	container *
shape:*
dtype0*
_output_shapes
:
Á
1layerfilter1PL_self_att_conv_head_0/biases/AssignAssign*layerfilter1PL_self_att_conv_head_0/biases<layerfilter1PL_self_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter1PL_self_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:
Ú
/layerfilter1PL_self_att_conv_head_0/biases/readIdentity*layerfilter1PL_self_att_conv_head_0/biases"/device:CPU:0*
_output_shapes
:*
T0*=
_class3
1/loc:@layerfilter1PL_self_att_conv_head_0/biases
Û
+layerfilter1PL_self_att_conv_head_0/BiasAddBiasAdd*layerfilter1PL_self_att_conv_head_0/Conv2D/layerfilter1PL_self_att_conv_head_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:
y
,layerfilter1PL_self_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

+layerfilter1PL_self_att_conv_head_0/bn/beta
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes
:*
	container 
¥
2layerfilter1PL_self_att_conv_head_0/bn/beta/AssignAssign+layerfilter1PL_self_att_conv_head_0/bn/beta,layerfilter1PL_self_att_conv_head_0/bn/Const*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
Î
0layerfilter1PL_self_att_conv_head_0/bn/beta/readIdentity+layerfilter1PL_self_att_conv_head_0/bn/beta*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/bn/beta*
_output_shapes
:*
T0
{
.layerfilter1PL_self_att_conv_head_0/bn/Const_1Const*
_output_shapes
:*
valueB*  ?*
dtype0

,layerfilter1PL_self_att_conv_head_0/bn/gamma
VariableV2*
_output_shapes
:*
	container *
shape:*
shared_name *
dtype0
ª
3layerfilter1PL_self_att_conv_head_0/bn/gamma/AssignAssign,layerfilter1PL_self_att_conv_head_0/bn/gamma.layerfilter1PL_self_att_conv_head_0/bn/Const_1*
use_locking(*
T0*?
_class5
31loc:@layerfilter1PL_self_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
Ñ
1layerfilter1PL_self_att_conv_head_0/bn/gamma/readIdentity,layerfilter1PL_self_att_conv_head_0/bn/gamma*
_output_shapes
:*
T0*?
_class5
31loc:@layerfilter1PL_self_att_conv_head_0/bn/gamma

Elayerfilter1PL_self_att_conv_head_0/bn/moments/mean/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0
ý
3layerfilter1PL_self_att_conv_head_0/bn/moments/meanMean+layerfilter1PL_self_att_conv_head_0/BiasAddElayerfilter1PL_self_att_conv_head_0/bn/moments/mean/reduction_indices*&
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
±
;layerfilter1PL_self_att_conv_head_0/bn/moments/StopGradientStopGradient3layerfilter1PL_self_att_conv_head_0/bn/moments/mean*&
_output_shapes
:*
T0
ð
@layerfilter1PL_self_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference+layerfilter1PL_self_att_conv_head_0/BiasAdd;layerfilter1PL_self_att_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:

Ilayerfilter1PL_self_att_conv_head_0/bn/moments/variance/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0

7layerfilter1PL_self_att_conv_head_0/bn/moments/varianceMean@layerfilter1PL_self_att_conv_head_0/bn/moments/SquaredDifferenceIlayerfilter1PL_self_att_conv_head_0/bn/moments/variance/reduction_indices*
	keep_dims(*

Tidx0*
T0*&
_output_shapes
:
´
6layerfilter1PL_self_att_conv_head_0/bn/moments/SqueezeSqueeze3layerfilter1PL_self_att_conv_head_0/bn/moments/mean*
squeeze_dims
 *
T0*
_output_shapes
:
º
8layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1Squeeze7layerfilter1PL_self_att_conv_head_0/bn/moments/variance*
_output_shapes
:*
squeeze_dims
 *
T0
}
2layerfilter1PL_self_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0


4layerfilter1PL_self_att_conv_head_0/bn/cond/switch_tIdentity4layerfilter1PL_self_att_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

4layerfilter1PL_self_att_conv_head_0/bn/cond/switch_fIdentity2layerfilter1PL_self_att_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
o
3layerfilter1PL_self_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
â
layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes
:
î
vlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shared_name *
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes
:
ä
}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
°
{layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ç
layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:
ó
xlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:*
shared_name *
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:
í
layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
·
}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
Jlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst5^layerfilter1PL_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ö
Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst5^layerfilter1PL_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¨
Xlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubZlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xJlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
Ò
Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subclayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1elayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes
:*
T0
Ô
alayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch{layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Ð
clayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch6layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze3layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id*I
_class?
=;loc:@layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::*
T0
º
Xlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulZlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Xlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
_output_shapes
:*
T0
Þ
Tlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub]layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Xlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:*
use_locking( *
T0
Ì
[layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage3layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Ø
\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst5^layerfilter1PL_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¬
Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xJlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
Ø
\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subelayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1glayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes
:*
T0
Û
clayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitch}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
elayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch8layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_13layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id*K
_classA
?=loc:@layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::*
T0
À
Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes
:*
T0
ç
Vlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub_layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:*
use_locking( *
T0
Ó
]layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage3layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
ü
Dlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpU^layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgW^layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
±
>layerfilter1PL_self_att_conv_head_0/bn/cond/control_dependencyIdentity4layerfilter1PL_self_att_conv_head_0/bn/cond/switch_tE^layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage*G
_class=
;9loc:@layerfilter1PL_self_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: *
T0

o
0layerfilter1PL_self_att_conv_head_0/bn/cond/NoOpNoOp5^layerfilter1PL_self_att_conv_head_0/bn/cond/switch_f

@layerfilter1PL_self_att_conv_head_0/bn/cond/control_dependency_1Identity4layerfilter1PL_self_att_conv_head_0/bn/cond/switch_f1^layerfilter1PL_self_att_conv_head_0/bn/cond/NoOp*
_output_shapes
: *
T0
*G
_class=
;9loc:@layerfilter1PL_self_att_conv_head_0/bn/cond/switch_f
è
1layerfilter1PL_self_att_conv_head_0/bn/cond/MergeMerge@layerfilter1PL_self_att_conv_head_0/bn/cond/control_dependency_1>layerfilter1PL_self_att_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 

4layerfilter1PL_self_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

6layerfilter1PL_self_att_conv_head_0/bn/cond_1/switch_tIdentity6layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch:1*
_output_shapes
: *
T0


6layerfilter1PL_self_att_conv_head_0/bn/cond_1/switch_fIdentity4layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
q
5layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

Ü
6layerfilter1PL_self_att_conv_head_0/bn/cond_1/IdentityIdentity?layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity/Switch:12^layerfilter1PL_self_att_conv_head_0/bn/cond/Merge*
_output_shapes
:*
T0
¬
=layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch6layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze5layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id*I
_class?
=;loc:@layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::*
T0
à
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity_1IdentityAlayerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:12^layerfilter1PL_self_att_conv_head_0/bn/cond/Merge*
_output_shapes
:*
T0
²
?layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch8layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_15layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id*
T0*K
_classA
?=loc:@layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
«
6layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_1Switch{layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read5layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::*
T0
°
6layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_2Switch}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read5layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::*
T0
Ü
3layerfilter1PL_self_att_conv_head_0/bn/cond_1/MergeMerge6layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_16layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
à
5layerfilter1PL_self_att_conv_head_0/bn/cond_1/Merge_1Merge6layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_28layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
{
6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
Ï
4layerfilter1PL_self_att_conv_head_0/bn/batchnorm/addAdd5layerfilter1PL_self_att_conv_head_0/bn/cond_1/Merge_16layerfilter1PL_self_att_conv_head_0/bn/batchnorm/add/y*
_output_shapes
:*
T0

6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/RsqrtRsqrt4layerfilter1PL_self_att_conv_head_0/bn/batchnorm/add*
_output_shapes
:*
T0
Ë
4layerfilter1PL_self_att_conv_head_0/bn/batchnorm/mulMul6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/Rsqrt1layerfilter1PL_self_att_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:
Ñ
6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/mul_1Mul+layerfilter1PL_self_att_conv_head_0/BiasAdd4layerfilter1PL_self_att_conv_head_0/bn/batchnorm/mul*
T0*&
_output_shapes
:
Í
6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/mul_2Mul3layerfilter1PL_self_att_conv_head_0/bn/cond_1/Merge4layerfilter1PL_self_att_conv_head_0/bn/batchnorm/mul*
_output_shapes
:*
T0
Ê
4layerfilter1PL_self_att_conv_head_0/bn/batchnorm/subSub0layerfilter1PL_self_att_conv_head_0/bn/beta/read6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:
Ü
6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/add_1Add6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/mul_14layerfilter1PL_self_att_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:

(layerfilter1PL_self_att_conv_head_0/ReluRelu6layerfilter1PL_self_att_conv_head_0/bn/batchnorm/add_1*&
_output_shapes
:*
T0
å
Llayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"      @      *>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
dtype0*
_output_shapes
:
Ï
Jlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *¾*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ï
Jlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *>*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
dtype0*
_output_shapes
: 
Ê
Tlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformLlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/shape*
dtype0*&
_output_shapes
:@*

seed *
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
seed2 
Ê
Jlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/subSubJlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/maxJlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
_output_shapes
: *
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights
ä
Jlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/mulMulTlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformJlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/sub*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*&
_output_shapes
:@
Ö
Flayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniformAddJlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/mulJlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*&
_output_shapes
:@
þ
+layerfilter1PL_neib_att_conv_head_0/weights
VariableV2"/device:CPU:0*
dtype0*&
_output_shapes
:@*
shared_name *>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
	container *
shape:@
Ú
2layerfilter1PL_neib_att_conv_head_0/weights/AssignAssign+layerfilter1PL_neib_att_conv_head_0/weightsFlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
validate_shape(*&
_output_shapes
:@
é
0layerfilter1PL_neib_att_conv_head_0/weights/readIdentity+layerfilter1PL_neib_att_conv_head_0/weights"/device:CPU:0*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*&
_output_shapes
:@

*layerfilter1PL_neib_att_conv_head_0/L2LossL2Loss0layerfilter1PL_neib_att_conv_head_0/weights/read*
T0*
_output_shapes
: 
v
1layerfilter1PL_neib_att_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
¶
/layerfilter1PL_neib_att_conv_head_0/weight_lossMul*layerfilter1PL_neib_att_conv_head_0/L2Loss1layerfilter1PL_neib_att_conv_head_0/weight_loss/y*
_output_shapes
: *
T0
¾
*layerfilter1PL_neib_att_conv_head_0/Conv2DConv2Dlayerfilter1PL_edgefea_0/Relu0layerfilter1PL_neib_att_conv_head_0/weights/read*
use_cudnn_on_gpu(*
explicit_paddings
 *
paddingVALID*&
_output_shapes
:
*
	dilations
*
T0*
strides
*
data_formatNHWC
È
<layerfilter1PL_neib_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *=
_class3
1/loc:@layerfilter1PL_neib_att_conv_head_0/biases*
dtype0*
_output_shapes
:
ä
*layerfilter1PL_neib_att_conv_head_0/biases
VariableV2"/device:CPU:0*
dtype0*
_output_shapes
:*
shared_name *=
_class3
1/loc:@layerfilter1PL_neib_att_conv_head_0/biases*
	container *
shape:
Á
1layerfilter1PL_neib_att_conv_head_0/biases/AssignAssign*layerfilter1PL_neib_att_conv_head_0/biases<layerfilter1PL_neib_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*=
_class3
1/loc:@layerfilter1PL_neib_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
Ú
/layerfilter1PL_neib_att_conv_head_0/biases/readIdentity*layerfilter1PL_neib_att_conv_head_0/biases"/device:CPU:0*=
_class3
1/loc:@layerfilter1PL_neib_att_conv_head_0/biases*
_output_shapes
:*
T0
Û
+layerfilter1PL_neib_att_conv_head_0/BiasAddBiasAdd*layerfilter1PL_neib_att_conv_head_0/Conv2D/layerfilter1PL_neib_att_conv_head_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:

y
,layerfilter1PL_neib_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

+layerfilter1PL_neib_att_conv_head_0/bn/beta
VariableV2*
_output_shapes
:*
	container *
shape:*
shared_name *
dtype0
¥
2layerfilter1PL_neib_att_conv_head_0/bn/beta/AssignAssign+layerfilter1PL_neib_att_conv_head_0/bn/beta,layerfilter1PL_neib_att_conv_head_0/bn/Const*
_output_shapes
:*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/bn/beta*
validate_shape(
Î
0layerfilter1PL_neib_att_conv_head_0/bn/beta/readIdentity+layerfilter1PL_neib_att_conv_head_0/bn/beta*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/bn/beta*
_output_shapes
:
{
.layerfilter1PL_neib_att_conv_head_0/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

,layerfilter1PL_neib_att_conv_head_0/bn/gamma
VariableV2*
_output_shapes
:*
	container *
shape:*
shared_name *
dtype0
ª
3layerfilter1PL_neib_att_conv_head_0/bn/gamma/AssignAssign,layerfilter1PL_neib_att_conv_head_0/bn/gamma.layerfilter1PL_neib_att_conv_head_0/bn/Const_1*
_output_shapes
:*
use_locking(*
T0*?
_class5
31loc:@layerfilter1PL_neib_att_conv_head_0/bn/gamma*
validate_shape(
Ñ
1layerfilter1PL_neib_att_conv_head_0/bn/gamma/readIdentity,layerfilter1PL_neib_att_conv_head_0/bn/gamma*
_output_shapes
:*
T0*?
_class5
31loc:@layerfilter1PL_neib_att_conv_head_0/bn/gamma

Elayerfilter1PL_neib_att_conv_head_0/bn/moments/mean/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0
ý
3layerfilter1PL_neib_att_conv_head_0/bn/moments/meanMean+layerfilter1PL_neib_att_conv_head_0/BiasAddElayerfilter1PL_neib_att_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*

Tidx0*
T0*&
_output_shapes
:
±
;layerfilter1PL_neib_att_conv_head_0/bn/moments/StopGradientStopGradient3layerfilter1PL_neib_att_conv_head_0/bn/moments/mean*&
_output_shapes
:*
T0
ð
@layerfilter1PL_neib_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference+layerfilter1PL_neib_att_conv_head_0/BiasAdd;layerfilter1PL_neib_att_conv_head_0/bn/moments/StopGradient*&
_output_shapes
:
*
T0

Ilayerfilter1PL_neib_att_conv_head_0/bn/moments/variance/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0

7layerfilter1PL_neib_att_conv_head_0/bn/moments/varianceMean@layerfilter1PL_neib_att_conv_head_0/bn/moments/SquaredDifferenceIlayerfilter1PL_neib_att_conv_head_0/bn/moments/variance/reduction_indices*&
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
´
6layerfilter1PL_neib_att_conv_head_0/bn/moments/SqueezeSqueeze3layerfilter1PL_neib_att_conv_head_0/bn/moments/mean*
T0*
_output_shapes
:*
squeeze_dims
 
º
8layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1Squeeze7layerfilter1PL_neib_att_conv_head_0/bn/moments/variance*
_output_shapes
:*
squeeze_dims
 *
T0
}
2layerfilter1PL_neib_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

4layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_tIdentity4layerfilter1PL_neib_att_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

4layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_fIdentity2layerfilter1PL_neib_att_conv_head_0/bn/cond/Switch*
_output_shapes
: *
T0

o
3layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
â
layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes
:
î
vlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes
:*
shared_name *
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:
ä
}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
°
{layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ç
layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:
ó
xlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
_output_shapes
:*
shared_name *
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0
í
layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
_output_shapes
:*
use_locking(*
T0*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
·
}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:*
T0*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
Æ
Jlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst5^layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ö
Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst5^layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¨
Xlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubZlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xJlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
Ò
Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subclayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1elayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_output_shapes
:
Ô
alayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch{layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Ð
clayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch6layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze3layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
º
Xlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulZlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Xlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes
:
Þ
Tlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub]layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Xlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
_output_shapes
:*
use_locking( *
T0*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
Ì
[layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage3layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Ø
\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst5^layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
¬
Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xJlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
Ø
\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subelayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1glayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes
:*
T0
Û
clayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitch}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::*
T0
Ö
elayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch8layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_13layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*K
_classA
?=loc:@layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
À
Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes
:*
T0
ç
Vlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub_layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
_output_shapes
:*
use_locking( *
T0*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
Ó
]layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage3layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
ü
Dlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpU^layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgW^layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
±
>layerfilter1PL_neib_att_conv_head_0/bn/cond/control_dependencyIdentity4layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_tE^layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage*G
_class=
;9loc:@layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: *
T0

o
0layerfilter1PL_neib_att_conv_head_0/bn/cond/NoOpNoOp5^layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_f

@layerfilter1PL_neib_att_conv_head_0/bn/cond/control_dependency_1Identity4layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_f1^layerfilter1PL_neib_att_conv_head_0/bn/cond/NoOp*G
_class=
;9loc:@layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_f*
_output_shapes
: *
T0

è
1layerfilter1PL_neib_att_conv_head_0/bn/cond/MergeMerge@layerfilter1PL_neib_att_conv_head_0/bn/cond/control_dependency_1>layerfilter1PL_neib_att_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 

4layerfilter1PL_neib_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 

6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/switch_tIdentity6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/switch_fIdentity4layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
q
5layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
Ü
6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/IdentityIdentity?layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity/Switch:12^layerfilter1PL_neib_att_conv_head_0/bn/cond/Merge*
_output_shapes
:*
T0
¬
=layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch6layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze5layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
à
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity_1IdentityAlayerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:12^layerfilter1PL_neib_att_conv_head_0/bn/cond/Merge*
_output_shapes
:*
T0
²
?layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch8layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_15layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id*K
_classA
?=loc:@layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::*
T0
«
6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_1Switch{layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read5layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::*
T0
°
6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_2Switch}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read5layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ü
3layerfilter1PL_neib_att_conv_head_0/bn/cond_1/MergeMerge6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_16layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
à
5layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Merge_1Merge6layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_28layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity_1*
N*
_output_shapes

:: *
T0
{
6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
Ï
4layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/addAdd5layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Merge_16layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/add/y*
_output_shapes
:*
T0

6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/RsqrtRsqrt4layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/add*
_output_shapes
:*
T0
Ë
4layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/mulMul6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/Rsqrt1layerfilter1PL_neib_att_conv_head_0/bn/gamma/read*
_output_shapes
:*
T0
Ñ
6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/mul_1Mul+layerfilter1PL_neib_att_conv_head_0/BiasAdd4layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/mul*&
_output_shapes
:
*
T0
Í
6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/mul_2Mul3layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Merge4layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:
Ê
4layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/subSub0layerfilter1PL_neib_att_conv_head_0/bn/beta/read6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/mul_2*
_output_shapes
:*
T0
Ü
6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/add_1Add6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/mul_14layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:


(layerfilter1PL_neib_att_conv_head_0/ReluRelu6layerfilter1PL_neib_att_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:


add_10Add(layerfilter1PL_self_att_conv_head_0/Relu(layerfilter1PL_neib_att_conv_head_0/Relu*&
_output_shapes
:
*
T0
i
transpose_7/permConst*
dtype0*
_output_shapes
:*%
valueB"             
p
transpose_7	Transposeadd_10transpose_7/perm*
T0*&
_output_shapes
:
*
Tperm0
f
LeakyRelu_1	LeakyRelutranspose_7*
T0*
alpha%ÍÌL>*&
_output_shapes
:

R
	Softmax_1SoftmaxLeakyRelu_1*&
_output_shapes
:
*
T0

MatMul_3BatchMatMulV2	Softmax_1layerfilter1PL_edgefea_0/Relu*&
_output_shapes
:@*
adj_x( *
adj_y( *
T0
®
/layerfilter1PL/BiasAdd/biases/Initializer/zerosConst*
valueB@*    *0
_class&
$"loc:@layerfilter1PL/BiasAdd/biases*
dtype0*
_output_shapes
:@
»
layerfilter1PL/BiasAdd/biases
VariableV2*
dtype0*
_output_shapes
:@*
shared_name *0
_class&
$"loc:@layerfilter1PL/BiasAdd/biases*
	container *
shape:@
þ
$layerfilter1PL/BiasAdd/biases/AssignAssignlayerfilter1PL/BiasAdd/biases/layerfilter1PL/BiasAdd/biases/Initializer/zeros*
_output_shapes
:@*
use_locking(*
T0*0
_class&
$"loc:@layerfilter1PL/BiasAdd/biases*
validate_shape(
¤
"layerfilter1PL/BiasAdd/biases/readIdentitylayerfilter1PL/BiasAdd/biases*
_output_shapes
:@*
T0*0
_class&
$"loc:@layerfilter1PL/BiasAdd/biases

layerfilter1PL/BiasAdd/BiasAddBiasAddMatMul_3"layerfilter1PL/BiasAdd/biases/read*&
_output_shapes
:@*
T0*
data_formatNHWC
_
Relu_1Relulayerfilter1PL/BiasAdd/BiasAdd*
T0*&
_output_shapes
:@
^
concat_3/concat_dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
M
concat_3IdentityRelu_1*&
_output_shapes
:@*
T0
\
ExpandDims_11/dimConst*
_output_shapes
: *
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0
x
ExpandDims_11
ExpandDimsPlaceholderExpandDims_11/dim*&
_output_shapes
:*

Tdim0*
T0
X
concat_4/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 

concat_4ConcatV2ExpandDims_11concat_3concat_4/axis*

Tidx0*
T0*
N*&
_output_shapes
:M
^
concat_5/concat_dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
d
concat_5Identitylayerfilter1PL_edgefea_0/Relu*
T0*&
_output_shapes
:
@
b
Max_2/reduction_indicesConst*
dtype0*
_output_shapes
: *
valueB :
þÿÿÿÿÿÿÿÿ
}
Max_2Maxconcat_5Max_2/reduction_indices*
T0*&
_output_shapes
:@*
	keep_dims(*

Tidx0
¯
1gapnet10/weights/Initializer/random_uniform/shapeConst*%
valueB"      M      *#
_class
loc:@gapnet10/weights*
dtype0*
_output_shapes
:

/gapnet10/weights/Initializer/random_uniform/minConst*
dtype0*
_output_shapes
: *
valueB
 *//¾*#
_class
loc:@gapnet10/weights

/gapnet10/weights/Initializer/random_uniform/maxConst*
valueB
 *//>*#
_class
loc:@gapnet10/weights*
dtype0*
_output_shapes
: 
ú
9gapnet10/weights/Initializer/random_uniform/RandomUniformRandomUniform1gapnet10/weights/Initializer/random_uniform/shape*
T0*#
_class
loc:@gapnet10/weights*
seed2 *
dtype0*'
_output_shapes
:M*

seed 
Þ
/gapnet10/weights/Initializer/random_uniform/subSub/gapnet10/weights/Initializer/random_uniform/max/gapnet10/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet10/weights*
_output_shapes
: 
ù
/gapnet10/weights/Initializer/random_uniform/mulMul9gapnet10/weights/Initializer/random_uniform/RandomUniform/gapnet10/weights/Initializer/random_uniform/sub*
T0*#
_class
loc:@gapnet10/weights*'
_output_shapes
:M
ë
+gapnet10/weights/Initializer/random_uniformAdd/gapnet10/weights/Initializer/random_uniform/mul/gapnet10/weights/Initializer/random_uniform/min*'
_output_shapes
:M*
T0*#
_class
loc:@gapnet10/weights
Ê
gapnet10/weights
VariableV2"/device:CPU:0*
dtype0*'
_output_shapes
:M*
shared_name *#
_class
loc:@gapnet10/weights*
	container *
shape:M
ï
gapnet10/weights/AssignAssigngapnet10/weights+gapnet10/weights/Initializer/random_uniform"/device:CPU:0*
T0*#
_class
loc:@gapnet10/weights*
validate_shape(*'
_output_shapes
:M*
use_locking(

gapnet10/weights/readIdentitygapnet10/weights"/device:CPU:0*
T0*#
_class
loc:@gapnet10/weights*'
_output_shapes
:M
Q
gapnet10/L2LossL2Lossgapnet10/weights/read*
T0*
_output_shapes
: 
[
gapnet10/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
e
gapnet10/weight_lossMulgapnet10/L2Lossgapnet10/weight_loss/y*
T0*
_output_shapes
: 
ô
gapnet10/Conv2DConv2Dconcat_4gapnet10/weights/read*
	dilations
*
T0*
data_formatNHWC*
strides
*
use_cudnn_on_gpu(*
explicit_paddings
 *
paddingVALID*'
_output_shapes
:

!gapnet10/biases/Initializer/ConstConst*
valueB*    *"
_class
loc:@gapnet10/biases*
dtype0*
_output_shapes	
:
°
gapnet10/biases
VariableV2"/device:CPU:0*"
_class
loc:@gapnet10/biases*
	container *
shape:*
dtype0*
_output_shapes	
:*
shared_name 
Ö
gapnet10/biases/AssignAssigngapnet10/biases!gapnet10/biases/Initializer/Const"/device:CPU:0*"
_class
loc:@gapnet10/biases*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

gapnet10/biases/readIdentitygapnet10/biases"/device:CPU:0*
T0*"
_class
loc:@gapnet10/biases*
_output_shapes	
:

gapnet10/BiasAddBiasAddgapnet10/Conv2Dgapnet10/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:
`
gapnet10/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
~
gapnet10/bn/beta
VariableV2*
_output_shapes	
:*
	container *
shape:*
shared_name *
dtype0
º
gapnet10/bn/beta/AssignAssigngapnet10/bn/betagapnet10/bn/Const*
use_locking(*
T0*#
_class
loc:@gapnet10/bn/beta*
validate_shape(*
_output_shapes	
:
~
gapnet10/bn/beta/readIdentitygapnet10/bn/beta*
T0*#
_class
loc:@gapnet10/bn/beta*
_output_shapes	
:
b
gapnet10/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:

gapnet10/bn/gamma
VariableV2*
_output_shapes	
:*
	container *
shape:*
shared_name *
dtype0
¿
gapnet10/bn/gamma/AssignAssigngapnet10/bn/gammagapnet10/bn/Const_1*
use_locking(*
T0*$
_class
loc:@gapnet10/bn/gamma*
validate_shape(*
_output_shapes	
:

gapnet10/bn/gamma/readIdentitygapnet10/bn/gamma*
T0*$
_class
loc:@gapnet10/bn/gamma*
_output_shapes	
:

*gapnet10/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
­
gapnet10/bn/moments/meanMeangapnet10/BiasAdd*gapnet10/bn/moments/mean/reduction_indices*'
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
|
 gapnet10/bn/moments/StopGradientStopGradientgapnet10/bn/moments/mean*
T0*'
_output_shapes
:
 
%gapnet10/bn/moments/SquaredDifferenceSquaredDifferencegapnet10/BiasAdd gapnet10/bn/moments/StopGradient*'
_output_shapes
:*
T0

.gapnet10/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ê
gapnet10/bn/moments/varianceMean%gapnet10/bn/moments/SquaredDifference.gapnet10/bn/moments/variance/reduction_indices*
T0*'
_output_shapes
:*
	keep_dims(*

Tidx0

gapnet10/bn/moments/SqueezeSqueezegapnet10/bn/moments/mean*
_output_shapes	
:*
squeeze_dims
 *
T0

gapnet10/bn/moments/Squeeze_1Squeezegapnet10/bn/moments/variance*
_output_shapes	
:*
squeeze_dims
 *
T0
b
gapnet10/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

a
gapnet10/bn/cond/switch_tIdentitygapnet10/bn/cond/Switch:1*
_output_shapes
: *
T0

_
gapnet10/bn/cond/switch_fIdentitygapnet10/bn/cond/Switch*
_output_shapes
: *
T0

T
gapnet10/bn/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

ö
Rgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes	
:*
valueB*    *S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0

@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
	container *
shape:*
dtype0*
_output_shapes	
:*
shared_name *S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage

Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageRgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
_output_shapes	
:*
use_locking(*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(

Egapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:*
T0
ú
Tgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes	
:

Bgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shared_name *U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes	
:

Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageTgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
_output_shapes	
:*
use_locking(*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(

Ggapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

/gapnet10/bn/cond/ExponentialMovingAverage/decayConst^gapnet10/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
 
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet10/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
×
=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x/gapnet10/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0

?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubHgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes	
:*
T0
³
Fgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchEgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet10/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
æ
Hgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet10/bn/moments/Squeezegapnet10/bn/cond/pred_id*"
_output_shapes
::*
T0*.
_class$
" loc:@gapnet10/bn/moments/Squeeze
ê
=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
_output_shapes	
:*
T0
×
9gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubBgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:*
use_locking( *
T0
«
@gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAveragegapnet10/bn/cond/pred_id*"
_output_shapes
::*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage
¢
Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet10/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Û
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubAgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x/gapnet10/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 

Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubJgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Lgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_output_shapes	
:
¹
Hgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchGgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet10/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
ì
Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchgapnet10/bn/moments/Squeeze_1gapnet10/bn/cond/pred_id*"
_output_shapes
::*
T0*0
_class&
$"loc:@gapnet10/bn/moments/Squeeze_1
ð
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulAgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_output_shapes	
:
ß
;gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubDgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
use_locking( *
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
±
Bgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet10/bn/cond/pred_id*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::*
T0
«
)gapnet10/bn/cond/ExponentialMovingAverageNoOp:^gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg<^gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
Å
#gapnet10/bn/cond/control_dependencyIdentitygapnet10/bn/cond/switch_t*^gapnet10/bn/cond/ExponentialMovingAverage*
T0
*,
_class"
 loc:@gapnet10/bn/cond/switch_t*
_output_shapes
: 
9
gapnet10/bn/cond/NoOpNoOp^gapnet10/bn/cond/switch_f
³
%gapnet10/bn/cond/control_dependency_1Identitygapnet10/bn/cond/switch_f^gapnet10/bn/cond/NoOp*,
_class"
 loc:@gapnet10/bn/cond/switch_f*
_output_shapes
: *
T0


gapnet10/bn/cond/MergeMerge%gapnet10/bn/cond/control_dependency_1#gapnet10/bn/cond/control_dependency*
_output_shapes
: : *
T0
*
N
d
gapnet10/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 
e
gapnet10/bn/cond_1/switch_tIdentitygapnet10/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
c
gapnet10/bn/cond_1/switch_fIdentitygapnet10/bn/cond_1/Switch*
_output_shapes
: *
T0

V
gapnet10/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 

gapnet10/bn/cond_1/IdentityIdentity$gapnet10/bn/cond_1/Identity/Switch:1^gapnet10/bn/cond/Merge*
_output_shapes	
:*
T0
Â
"gapnet10/bn/cond_1/Identity/SwitchSwitchgapnet10/bn/moments/Squeezegapnet10/bn/cond_1/pred_id*.
_class$
" loc:@gapnet10/bn/moments/Squeeze*"
_output_shapes
::*
T0

gapnet10/bn/cond_1/Identity_1Identity&gapnet10/bn/cond_1/Identity_1/Switch:1^gapnet10/bn/cond/Merge*
_output_shapes	
:*
T0
È
$gapnet10/bn/cond_1/Identity_1/SwitchSwitchgapnet10/bn/moments/Squeeze_1gapnet10/bn/cond_1/pred_id*"
_output_shapes
::*
T0*0
_class&
$"loc:@gapnet10/bn/moments/Squeeze_1

gapnet10/bn/cond_1/Switch_1SwitchEgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet10/bn/cond_1/pred_id*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

gapnet10/bn/cond_1/Switch_2SwitchGgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet10/bn/cond_1/pred_id*"
_output_shapes
::*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage

gapnet10/bn/cond_1/MergeMergegapnet10/bn/cond_1/Switch_1gapnet10/bn/cond_1/Identity*
T0*
N*
_output_shapes
	:: 

gapnet10/bn/cond_1/Merge_1Mergegapnet10/bn/cond_1/Switch_2gapnet10/bn/cond_1/Identity_1*
N*
_output_shapes
	:: *
T0
`
gapnet10/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 

gapnet10/bn/batchnorm/addAddgapnet10/bn/cond_1/Merge_1gapnet10/bn/batchnorm/add/y*
_output_shapes	
:*
T0
e
gapnet10/bn/batchnorm/RsqrtRsqrtgapnet10/bn/batchnorm/add*
_output_shapes	
:*
T0
{
gapnet10/bn/batchnorm/mulMulgapnet10/bn/batchnorm/Rsqrtgapnet10/bn/gamma/read*
T0*
_output_shapes	
:

gapnet10/bn/batchnorm/mul_1Mulgapnet10/BiasAddgapnet10/bn/batchnorm/mul*'
_output_shapes
:*
T0
}
gapnet10/bn/batchnorm/mul_2Mulgapnet10/bn/cond_1/Mergegapnet10/bn/batchnorm/mul*
_output_shapes	
:*
T0
z
gapnet10/bn/batchnorm/subSubgapnet10/bn/beta/readgapnet10/bn/batchnorm/mul_2*
_output_shapes	
:*
T0

gapnet10/bn/batchnorm/add_1Addgapnet10/bn/batchnorm/mul_1gapnet10/bn/batchnorm/sub*'
_output_shapes
:*
T0
d
gapnet10/ReluRelugapnet10/bn/batchnorm/add_1*
T0*'
_output_shapes
:
³
3gapnet11PL/weights/Initializer/random_uniform/shapeConst*%
valueB"            *%
_class
loc:@gapnet11PL/weights*
dtype0*
_output_shapes
:

1gapnet11PL/weights/Initializer/random_uniform/minConst*
valueB
 *qÄ¾*%
_class
loc:@gapnet11PL/weights*
dtype0*
_output_shapes
: 

1gapnet11PL/weights/Initializer/random_uniform/maxConst*
valueB
 *qÄ>*%
_class
loc:@gapnet11PL/weights*
dtype0*
_output_shapes
: 

;gapnet11PL/weights/Initializer/random_uniform/RandomUniformRandomUniform3gapnet11PL/weights/Initializer/random_uniform/shape*
dtype0*(
_output_shapes
:*

seed *
T0*%
_class
loc:@gapnet11PL/weights*
seed2 
æ
1gapnet11PL/weights/Initializer/random_uniform/subSub1gapnet11PL/weights/Initializer/random_uniform/max1gapnet11PL/weights/Initializer/random_uniform/min*%
_class
loc:@gapnet11PL/weights*
_output_shapes
: *
T0

1gapnet11PL/weights/Initializer/random_uniform/mulMul;gapnet11PL/weights/Initializer/random_uniform/RandomUniform1gapnet11PL/weights/Initializer/random_uniform/sub*(
_output_shapes
:*
T0*%
_class
loc:@gapnet11PL/weights
ô
-gapnet11PL/weights/Initializer/random_uniformAdd1gapnet11PL/weights/Initializer/random_uniform/mul1gapnet11PL/weights/Initializer/random_uniform/min*%
_class
loc:@gapnet11PL/weights*(
_output_shapes
:*
T0
Ð
gapnet11PL/weights
VariableV2"/device:CPU:0*
shape:*
dtype0*(
_output_shapes
:*
shared_name *%
_class
loc:@gapnet11PL/weights*
	container 
ø
gapnet11PL/weights/AssignAssigngapnet11PL/weights-gapnet11PL/weights/Initializer/random_uniform"/device:CPU:0*
use_locking(*
T0*%
_class
loc:@gapnet11PL/weights*
validate_shape(*(
_output_shapes
:
 
gapnet11PL/weights/readIdentitygapnet11PL/weights"/device:CPU:0*
T0*%
_class
loc:@gapnet11PL/weights*(
_output_shapes
:
U
gapnet11PL/L2LossL2Lossgapnet11PL/weights/read*
_output_shapes
: *
T0
]
gapnet11PL/weight_loss/yConst*
_output_shapes
: *
valueB
 *    *
dtype0
k
gapnet11PL/weight_lossMulgapnet11PL/L2Lossgapnet11PL/weight_loss/y*
_output_shapes
: *
T0
ý
gapnet11PL/Conv2DConv2Dgapnet10/Relugapnet11PL/weights/read*
paddingVALID*'
_output_shapes
:*
	dilations
*
T0*
data_formatNHWC*
strides
*
use_cudnn_on_gpu(*
explicit_paddings
 

#gapnet11PL/biases/Initializer/ConstConst*
valueB*    *$
_class
loc:@gapnet11PL/biases*
dtype0*
_output_shapes	
:
´
gapnet11PL/biases
VariableV2"/device:CPU:0*
	container *
shape:*
dtype0*
_output_shapes	
:*
shared_name *$
_class
loc:@gapnet11PL/biases
Þ
gapnet11PL/biases/AssignAssigngapnet11PL/biases#gapnet11PL/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*$
_class
loc:@gapnet11PL/biases*
validate_shape(*
_output_shapes	
:

gapnet11PL/biases/readIdentitygapnet11PL/biases"/device:CPU:0*$
_class
loc:@gapnet11PL/biases*
_output_shapes	
:*
T0

gapnet11PL/BiasAddBiasAddgapnet11PL/Conv2Dgapnet11PL/biases/read*'
_output_shapes
:*
T0*
data_formatNHWC
b
gapnet11PL/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:

gapnet11PL/bn/beta
VariableV2*
shape:*
shared_name *
dtype0*
_output_shapes	
:*
	container 
Â
gapnet11PL/bn/beta/AssignAssigngapnet11PL/bn/betagapnet11PL/bn/Const*
use_locking(*
T0*%
_class
loc:@gapnet11PL/bn/beta*
validate_shape(*
_output_shapes	
:

gapnet11PL/bn/beta/readIdentitygapnet11PL/bn/beta*
T0*%
_class
loc:@gapnet11PL/bn/beta*
_output_shapes	
:
d
gapnet11PL/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:

gapnet11PL/bn/gamma
VariableV2*
_output_shapes	
:*
	container *
shape:*
shared_name *
dtype0
Ç
gapnet11PL/bn/gamma/AssignAssigngapnet11PL/bn/gammagapnet11PL/bn/Const_1*
use_locking(*
T0*&
_class
loc:@gapnet11PL/bn/gamma*
validate_shape(*
_output_shapes	
:

gapnet11PL/bn/gamma/readIdentitygapnet11PL/bn/gamma*&
_class
loc:@gapnet11PL/bn/gamma*
_output_shapes	
:*
T0

,gapnet11PL/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
³
gapnet11PL/bn/moments/meanMeangapnet11PL/BiasAdd,gapnet11PL/bn/moments/mean/reduction_indices*'
_output_shapes
:*
	keep_dims(*

Tidx0*
T0

"gapnet11PL/bn/moments/StopGradientStopGradientgapnet11PL/bn/moments/mean*'
_output_shapes
:*
T0
¦
'gapnet11PL/bn/moments/SquaredDifferenceSquaredDifferencegapnet11PL/BiasAdd"gapnet11PL/bn/moments/StopGradient*
T0*'
_output_shapes
:

0gapnet11PL/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ð
gapnet11PL/bn/moments/varianceMean'gapnet11PL/bn/moments/SquaredDifference0gapnet11PL/bn/moments/variance/reduction_indices*'
_output_shapes
:*
	keep_dims(*

Tidx0*
T0

gapnet11PL/bn/moments/SqueezeSqueezegapnet11PL/bn/moments/mean*
squeeze_dims
 *
T0*
_output_shapes	
:

gapnet11PL/bn/moments/Squeeze_1Squeezegapnet11PL/bn/moments/variance*
squeeze_dims
 *
T0*
_output_shapes	
:
d
gapnet11PL/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 
e
gapnet11PL/bn/cond/switch_tIdentitygapnet11PL/bn/cond/Switch:1*
T0
*
_output_shapes
: 
c
gapnet11PL/bn/cond/switch_fIdentitygapnet11PL/bn/cond/Switch*
_output_shapes
: *
T0

V
gapnet11PL/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
þ
Vgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes	
:*
valueB*    *W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0

Dgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
_output_shapes	
:*
shared_name *W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:*
dtype0

Kgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverageVgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

Igapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
T0*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Xgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes	
:

Fgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes	
:*
shared_name *Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:
£
Mgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverageXgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
_output_shapes	
:*
use_locking(*
T0*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
 
Kgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:*
T0

1gapnet11PL/bn/cond/ExponentialMovingAverage/decayConst^gapnet11PL/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
¤
Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet11PL/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
Ý
?gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubAgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x1gapnet11PL/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 

Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubJgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Lgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_output_shapes	
:
¿
Hgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchIgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet11PL/bn/cond/pred_id*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::*
T0
î
Jgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet11PL/bn/moments/Squeezegapnet11PL/bn/cond/pred_id*
T0*0
_class&
$"loc:@gapnet11PL/bn/moments/Squeeze*"
_output_shapes
::
ð
?gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulAgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1?gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
_output_shapes	
:*
T0
á
;gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubDgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1?gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
use_locking( *
T0*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
·
Bgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAveragegapnet11PL/bn/cond/pred_id*
T0*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
¦
Cgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet11PL/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
á
Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubCgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x1gapnet11PL/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0

Cgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubLgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Ngapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_output_shapes	
:
Å
Jgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchKgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet11PL/bn/cond/pred_id*
T0*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
ô
Lgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchgapnet11PL/bn/moments/Squeeze_1gapnet11PL/bn/cond/pred_id*
T0*2
_class(
&$loc:@gapnet11PL/bn/moments/Squeeze_1*"
_output_shapes
::
ö
Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulCgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_output_shapes	
:
é
=gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubFgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
use_locking( *
T0*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
½
Dgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet11PL/bn/cond/pred_id*"
_output_shapes
::*
T0*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage
±
+gapnet11PL/bn/cond/ExponentialMovingAverageNoOp<^gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg>^gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
Í
%gapnet11PL/bn/cond/control_dependencyIdentitygapnet11PL/bn/cond/switch_t,^gapnet11PL/bn/cond/ExponentialMovingAverage*
T0
*.
_class$
" loc:@gapnet11PL/bn/cond/switch_t*
_output_shapes
: 
=
gapnet11PL/bn/cond/NoOpNoOp^gapnet11PL/bn/cond/switch_f
»
'gapnet11PL/bn/cond/control_dependency_1Identitygapnet11PL/bn/cond/switch_f^gapnet11PL/bn/cond/NoOp*
T0
*.
_class$
" loc:@gapnet11PL/bn/cond/switch_f*
_output_shapes
: 

gapnet11PL/bn/cond/MergeMerge'gapnet11PL/bn/cond/control_dependency_1%gapnet11PL/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
f
gapnet11PL/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

i
gapnet11PL/bn/cond_1/switch_tIdentitygapnet11PL/bn/cond_1/Switch:1*
_output_shapes
: *
T0

g
gapnet11PL/bn/cond_1/switch_fIdentitygapnet11PL/bn/cond_1/Switch*
T0
*
_output_shapes
: 
X
gapnet11PL/bn/cond_1/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0


gapnet11PL/bn/cond_1/IdentityIdentity&gapnet11PL/bn/cond_1/Identity/Switch:1^gapnet11PL/bn/cond/Merge*
_output_shapes	
:*
T0
Ê
$gapnet11PL/bn/cond_1/Identity/SwitchSwitchgapnet11PL/bn/moments/Squeezegapnet11PL/bn/cond_1/pred_id*0
_class&
$"loc:@gapnet11PL/bn/moments/Squeeze*"
_output_shapes
::*
T0

gapnet11PL/bn/cond_1/Identity_1Identity(gapnet11PL/bn/cond_1/Identity_1/Switch:1^gapnet11PL/bn/cond/Merge*
_output_shapes	
:*
T0
Ð
&gapnet11PL/bn/cond_1/Identity_1/SwitchSwitchgapnet11PL/bn/moments/Squeeze_1gapnet11PL/bn/cond_1/pred_id*"
_output_shapes
::*
T0*2
_class(
&$loc:@gapnet11PL/bn/moments/Squeeze_1

gapnet11PL/bn/cond_1/Switch_1SwitchIgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet11PL/bn/cond_1/pred_id*
T0*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

gapnet11PL/bn/cond_1/Switch_2SwitchKgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet11PL/bn/cond_1/pred_id*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::*
T0

gapnet11PL/bn/cond_1/MergeMergegapnet11PL/bn/cond_1/Switch_1gapnet11PL/bn/cond_1/Identity*
N*
_output_shapes
	:: *
T0

gapnet11PL/bn/cond_1/Merge_1Mergegapnet11PL/bn/cond_1/Switch_2gapnet11PL/bn/cond_1/Identity_1*
_output_shapes
	:: *
T0*
N
b
gapnet11PL/bn/batchnorm/add/yConst*
_output_shapes
: *
valueB
 *o:*
dtype0

gapnet11PL/bn/batchnorm/addAddgapnet11PL/bn/cond_1/Merge_1gapnet11PL/bn/batchnorm/add/y*
T0*
_output_shapes	
:
i
gapnet11PL/bn/batchnorm/RsqrtRsqrtgapnet11PL/bn/batchnorm/add*
T0*
_output_shapes	
:

gapnet11PL/bn/batchnorm/mulMulgapnet11PL/bn/batchnorm/Rsqrtgapnet11PL/bn/gamma/read*
T0*
_output_shapes	
:

gapnet11PL/bn/batchnorm/mul_1Mulgapnet11PL/BiasAddgapnet11PL/bn/batchnorm/mul*
T0*'
_output_shapes
:

gapnet11PL/bn/batchnorm/mul_2Mulgapnet11PL/bn/cond_1/Mergegapnet11PL/bn/batchnorm/mul*
T0*
_output_shapes	
:

gapnet11PL/bn/batchnorm/subSubgapnet11PL/bn/beta/readgapnet11PL/bn/batchnorm/mul_2*
_output_shapes	
:*
T0

gapnet11PL/bn/batchnorm/add_1Addgapnet11PL/bn/batchnorm/mul_1gapnet11PL/bn/batchnorm/sub*'
_output_shapes
:*
T0
h
gapnet11PL/ReluRelugapnet11PL/bn/batchnorm/add_1*'
_output_shapes
:*
T0
b
Max_3/reduction_indicesConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 

Max_3Maxgapnet11PL/ReluMax_3/reduction_indices*
T0*'
_output_shapes
:*
	keep_dims(*

Tidx0
X
concat_6/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
¢
concat_6ConcatV2gapnet00/ReluMax_1gapnet10/ReluMax_3MaxMax_2concat_6/axis*
T0*
N*'
_output_shapes
:*

Tidx0
©
.aggPL/weights/Initializer/random_uniform/shapeConst*%
valueB"           * 
_class
loc:@aggPL/weights*
dtype0*
_output_shapes
:

,aggPL/weights/Initializer/random_uniform/minConst*
valueB
 *¥)³½* 
_class
loc:@aggPL/weights*
dtype0*
_output_shapes
: 

,aggPL/weights/Initializer/random_uniform/maxConst*
valueB
 *¥)³=* 
_class
loc:@aggPL/weights*
dtype0*
_output_shapes
: 
ò
6aggPL/weights/Initializer/random_uniform/RandomUniformRandomUniform.aggPL/weights/Initializer/random_uniform/shape*
dtype0*(
_output_shapes
:*

seed *
T0* 
_class
loc:@aggPL/weights*
seed2 
Ò
,aggPL/weights/Initializer/random_uniform/subSub,aggPL/weights/Initializer/random_uniform/max,aggPL/weights/Initializer/random_uniform/min*
_output_shapes
: *
T0* 
_class
loc:@aggPL/weights
î
,aggPL/weights/Initializer/random_uniform/mulMul6aggPL/weights/Initializer/random_uniform/RandomUniform,aggPL/weights/Initializer/random_uniform/sub* 
_class
loc:@aggPL/weights*(
_output_shapes
:*
T0
à
(aggPL/weights/Initializer/random_uniformAdd,aggPL/weights/Initializer/random_uniform/mul,aggPL/weights/Initializer/random_uniform/min* 
_class
loc:@aggPL/weights*(
_output_shapes
:*
T0
Æ
aggPL/weights
VariableV2"/device:CPU:0* 
_class
loc:@aggPL/weights*
	container *
shape:*
dtype0*(
_output_shapes
:*
shared_name 
ä
aggPL/weights/AssignAssignaggPL/weights(aggPL/weights/Initializer/random_uniform"/device:CPU:0*(
_output_shapes
:*
use_locking(*
T0* 
_class
loc:@aggPL/weights*
validate_shape(

aggPL/weights/readIdentityaggPL/weights"/device:CPU:0* 
_class
loc:@aggPL/weights*(
_output_shapes
:*
T0
K
aggPL/L2LossL2LossaggPL/weights/read*
T0*
_output_shapes
: 
X
aggPL/weight_loss/yConst*
_output_shapes
: *
valueB
 *    *
dtype0
\
aggPL/weight_lossMulaggPL/L2LossaggPL/weight_loss/y*
_output_shapes
: *
T0
î
aggPL/Conv2DConv2Dconcat_6aggPL/weights/read*
paddingVALID*'
_output_shapes
:*
	dilations
*
T0*
strides
*
data_formatNHWC*
use_cudnn_on_gpu(*
explicit_paddings
 

aggPL/biases/Initializer/ConstConst*
_output_shapes	
:*
valueB*    *
_class
loc:@aggPL/biases*
dtype0
ª
aggPL/biases
VariableV2"/device:CPU:0*
shared_name *
_class
loc:@aggPL/biases*
	container *
shape:*
dtype0*
_output_shapes	
:
Ê
aggPL/biases/AssignAssignaggPL/biasesaggPL/biases/Initializer/Const"/device:CPU:0*
_output_shapes	
:*
use_locking(*
T0*
_class
loc:@aggPL/biases*
validate_shape(

aggPL/biases/readIdentityaggPL/biases"/device:CPU:0*
_class
loc:@aggPL/biases*
_output_shapes	
:*
T0

aggPL/BiasAddBiasAddaggPL/Conv2DaggPL/biases/read*'
_output_shapes
:*
T0*
data_formatNHWC
]
aggPL/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
{
aggPL/bn/beta
VariableV2*
shared_name *
dtype0*
_output_shapes	
:*
	container *
shape:
®
aggPL/bn/beta/AssignAssignaggPL/bn/betaaggPL/bn/Const*
_output_shapes	
:*
use_locking(*
T0* 
_class
loc:@aggPL/bn/beta*
validate_shape(
u
aggPL/bn/beta/readIdentityaggPL/bn/beta*
T0* 
_class
loc:@aggPL/bn/beta*
_output_shapes	
:
_
aggPL/bn/Const_1Const*
_output_shapes	
:*
valueB*  ?*
dtype0
|
aggPL/bn/gamma
VariableV2*
dtype0*
_output_shapes	
:*
	container *
shape:*
shared_name 
³
aggPL/bn/gamma/AssignAssignaggPL/bn/gammaaggPL/bn/Const_1*!
_class
loc:@aggPL/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
x
aggPL/bn/gamma/readIdentityaggPL/bn/gamma*
T0*!
_class
loc:@aggPL/bn/gamma*
_output_shapes	
:
|
'aggPL/bn/moments/mean/reduction_indicesConst*
_output_shapes
:*!
valueB"          *
dtype0
¤
aggPL/bn/moments/meanMeanaggPL/BiasAdd'aggPL/bn/moments/mean/reduction_indices*'
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
v
aggPL/bn/moments/StopGradientStopGradientaggPL/bn/moments/mean*
T0*'
_output_shapes
:

"aggPL/bn/moments/SquaredDifferenceSquaredDifferenceaggPL/BiasAddaggPL/bn/moments/StopGradient*
T0*'
_output_shapes
:

+aggPL/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Á
aggPL/bn/moments/varianceMean"aggPL/bn/moments/SquaredDifference+aggPL/bn/moments/variance/reduction_indices*'
_output_shapes
:*
	keep_dims(*

Tidx0*
T0
y
aggPL/bn/moments/SqueezeSqueezeaggPL/bn/moments/mean*
T0*
_output_shapes	
:*
squeeze_dims
 

aggPL/bn/moments/Squeeze_1SqueezeaggPL/bn/moments/variance*
T0*
_output_shapes	
:*
squeeze_dims
 
_
aggPL/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

[
aggPL/bn/cond/switch_tIdentityaggPL/bn/cond/Switch:1*
_output_shapes
: *
T0

Y
aggPL/bn/cond/switch_fIdentityaggPL/bn/cond/Switch*
_output_shapes
: *
T0

Q
aggPL/bn/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

ê
LaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes	
:
÷
:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes	
:*
shared_name *M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:
ó
AaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverageLaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
_output_shapes	
:*
use_locking(*
T0*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(
ü
?aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
T0*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
î
NaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes	
:*
valueB*    *O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0
û
<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes	
:*
shared_name 
û
CaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssign<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageNaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

AaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentity<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:*
T0*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage

,aggPL/bn/cond/ExponentialMovingAverage/decayConst^aggPL/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 

<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^aggPL/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Î
:aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x,aggPL/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
ù
<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubEaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1GaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes	
:*
T0
¡
CaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch?aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/readaggPL/bn/cond/pred_id*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::*
T0
Ú
EaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1SwitchaggPL/bn/moments/SqueezeaggPL/bn/cond/pred_id*+
_class!
loc:@aggPL/bn/moments/Squeeze*"
_output_shapes
::*
T0
á
:aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes	
:
È
6aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub?aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1:aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:*
use_locking( *
T0

=aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverageaggPL/bn/cond/pred_id*
T0*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^aggPL/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
Ò
<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x,aggPL/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
ÿ
>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubGaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1IaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_output_shapes	
:
§
EaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchAaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/readaggPL/bn/cond/pred_id*"
_output_shapes
::*
T0*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage
à
GaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1SwitchaggPL/bn/moments/Squeeze_1aggPL/bn/cond/pred_id*-
_class#
!loc:@aggPL/bn/moments/Squeeze_1*"
_output_shapes
::*
T0
ç
<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_output_shapes	
:
Ð
8aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubAaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:*
use_locking( *
T0

?aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitch<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageaggPL/bn/cond/pred_id*
T0*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
¢
&aggPL/bn/cond/ExponentialMovingAverageNoOp7^aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg9^aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
¹
 aggPL/bn/cond/control_dependencyIdentityaggPL/bn/cond/switch_t'^aggPL/bn/cond/ExponentialMovingAverage*
T0
*)
_class
loc:@aggPL/bn/cond/switch_t*
_output_shapes
: 
3
aggPL/bn/cond/NoOpNoOp^aggPL/bn/cond/switch_f
§
"aggPL/bn/cond/control_dependency_1IdentityaggPL/bn/cond/switch_f^aggPL/bn/cond/NoOp*
T0
*)
_class
loc:@aggPL/bn/cond/switch_f*
_output_shapes
: 

aggPL/bn/cond/MergeMerge"aggPL/bn/cond/control_dependency_1 aggPL/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
a
aggPL/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 
_
aggPL/bn/cond_1/switch_tIdentityaggPL/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
]
aggPL/bn/cond_1/switch_fIdentityaggPL/bn/cond_1/Switch*
_output_shapes
: *
T0

S
aggPL/bn/cond_1/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0


aggPL/bn/cond_1/IdentityIdentity!aggPL/bn/cond_1/Identity/Switch:1^aggPL/bn/cond/Merge*
_output_shapes	
:*
T0
¶
aggPL/bn/cond_1/Identity/SwitchSwitchaggPL/bn/moments/SqueezeaggPL/bn/cond_1/pred_id*+
_class!
loc:@aggPL/bn/moments/Squeeze*"
_output_shapes
::*
T0

aggPL/bn/cond_1/Identity_1Identity#aggPL/bn/cond_1/Identity_1/Switch:1^aggPL/bn/cond/Merge*
_output_shapes	
:*
T0
¼
!aggPL/bn/cond_1/Identity_1/SwitchSwitchaggPL/bn/moments/Squeeze_1aggPL/bn/cond_1/pred_id*"
_output_shapes
::*
T0*-
_class#
!loc:@aggPL/bn/moments/Squeeze_1
ø
aggPL/bn/cond_1/Switch_1Switch?aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/readaggPL/bn/cond_1/pred_id*
T0*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
ü
aggPL/bn/cond_1/Switch_2SwitchAaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/readaggPL/bn/cond_1/pred_id*
T0*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::

aggPL/bn/cond_1/MergeMergeaggPL/bn/cond_1/Switch_1aggPL/bn/cond_1/Identity*
N*
_output_shapes
	:: *
T0

aggPL/bn/cond_1/Merge_1MergeaggPL/bn/cond_1/Switch_2aggPL/bn/cond_1/Identity_1*
T0*
N*
_output_shapes
	:: 
]
aggPL/bn/batchnorm/add/yConst*
_output_shapes
: *
valueB
 *o:*
dtype0
v
aggPL/bn/batchnorm/addAddaggPL/bn/cond_1/Merge_1aggPL/bn/batchnorm/add/y*
_output_shapes	
:*
T0
_
aggPL/bn/batchnorm/RsqrtRsqrtaggPL/bn/batchnorm/add*
_output_shapes	
:*
T0
r
aggPL/bn/batchnorm/mulMulaggPL/bn/batchnorm/RsqrtaggPL/bn/gamma/read*
T0*
_output_shapes	
:
x
aggPL/bn/batchnorm/mul_1MulaggPL/BiasAddaggPL/bn/batchnorm/mul*
T0*'
_output_shapes
:
t
aggPL/bn/batchnorm/mul_2MulaggPL/bn/cond_1/MergeaggPL/bn/batchnorm/mul*
_output_shapes	
:*
T0
q
aggPL/bn/batchnorm/subSubaggPL/bn/beta/readaggPL/bn/batchnorm/mul_2*
T0*
_output_shapes	
:

aggPL/bn/batchnorm/add_1AddaggPL/bn/batchnorm/mul_1aggPL/bn/batchnorm/sub*
T0*'
_output_shapes
:
^

aggPL/ReluReluaggPL/bn/batchnorm/add_1*
T0*'
_output_shapes
:
¯
avgpoolPL/avgpoolPLMaxPool
aggPL/Relu*
T0*
data_formatNHWC*
strides
*
ksize
*
paddingVALID*'
_output_shapes
:
`
Reshape_4/shapeConst*
valueB"   ÿÿÿÿ*
dtype0*
_output_shapes
:
r
	Reshape_4ReshapeavgpoolPL/avgpoolPLReshape_4/shape*
Tshape0*
_output_shapes
:	*
T0
¡
.PLfc1/weights/Initializer/random_uniform/shapeConst*
valueB"      * 
_class
loc:@PLfc1/weights*
dtype0*
_output_shapes
:

,PLfc1/weights/Initializer/random_uniform/minConst*
valueB
 *   ¾* 
_class
loc:@PLfc1/weights*
dtype0*
_output_shapes
: 

,PLfc1/weights/Initializer/random_uniform/maxConst*
_output_shapes
: *
valueB
 *   >* 
_class
loc:@PLfc1/weights*
dtype0
ê
6PLfc1/weights/Initializer/random_uniform/RandomUniformRandomUniform.PLfc1/weights/Initializer/random_uniform/shape*

seed *
T0* 
_class
loc:@PLfc1/weights*
seed2 *
dtype0* 
_output_shapes
:

Ò
,PLfc1/weights/Initializer/random_uniform/subSub,PLfc1/weights/Initializer/random_uniform/max,PLfc1/weights/Initializer/random_uniform/min* 
_class
loc:@PLfc1/weights*
_output_shapes
: *
T0
æ
,PLfc1/weights/Initializer/random_uniform/mulMul6PLfc1/weights/Initializer/random_uniform/RandomUniform,PLfc1/weights/Initializer/random_uniform/sub* 
_class
loc:@PLfc1/weights* 
_output_shapes
:
*
T0
Ø
(PLfc1/weights/Initializer/random_uniformAdd,PLfc1/weights/Initializer/random_uniform/mul,PLfc1/weights/Initializer/random_uniform/min*
T0* 
_class
loc:@PLfc1/weights* 
_output_shapes
:

¶
PLfc1/weights
VariableV2"/device:CPU:0*
shape:
*
dtype0* 
_output_shapes
:
*
shared_name * 
_class
loc:@PLfc1/weights*
	container 
Ü
PLfc1/weights/AssignAssignPLfc1/weights(PLfc1/weights/Initializer/random_uniform"/device:CPU:0*
use_locking(*
T0* 
_class
loc:@PLfc1/weights*
validate_shape(* 
_output_shapes
:


PLfc1/weights/readIdentityPLfc1/weights"/device:CPU:0* 
_output_shapes
:
*
T0* 
_class
loc:@PLfc1/weights
K
PLfc1/L2LossL2LossPLfc1/weights/read*
_output_shapes
: *
T0
X
PLfc1/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
\
PLfc1/weight_lossMulPLfc1/L2LossPLfc1/weight_loss/y*
_output_shapes
: *
T0

PLfc1/MatMulMatMul	Reshape_4PLfc1/weights/read*
_output_shapes
:	*
transpose_a( *
transpose_b( *
T0

PLfc1/biases/Initializer/ConstConst*
valueB*    *
_class
loc:@PLfc1/biases*
dtype0*
_output_shapes	
:
ª
PLfc1/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
_output_shapes	
:*
shared_name *
_class
loc:@PLfc1/biases*
	container 
Ê
PLfc1/biases/AssignAssignPLfc1/biasesPLfc1/biases/Initializer/Const"/device:CPU:0*
_output_shapes	
:*
use_locking(*
T0*
_class
loc:@PLfc1/biases*
validate_shape(

PLfc1/biases/readIdentityPLfc1/biases"/device:CPU:0*
_output_shapes	
:*
T0*
_class
loc:@PLfc1/biases
z
PLfc1/BiasAddBiasAddPLfc1/MatMulPLfc1/biases/read*
data_formatNHWC*
_output_shapes
:	*
T0
]
PLfc1/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
{
PLfc1/bn/beta
VariableV2*
_output_shapes	
:*
	container *
shape:*
shared_name *
dtype0
®
PLfc1/bn/beta/AssignAssignPLfc1/bn/betaPLfc1/bn/Const* 
_class
loc:@PLfc1/bn/beta*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
u
PLfc1/bn/beta/readIdentityPLfc1/bn/beta* 
_class
loc:@PLfc1/bn/beta*
_output_shapes	
:*
T0
_
PLfc1/bn/Const_1Const*
_output_shapes	
:*
valueB*  ?*
dtype0
|
PLfc1/bn/gamma
VariableV2*
shared_name *
dtype0*
_output_shapes	
:*
	container *
shape:
³
PLfc1/bn/gamma/AssignAssignPLfc1/bn/gammaPLfc1/bn/Const_1*
use_locking(*
T0*!
_class
loc:@PLfc1/bn/gamma*
validate_shape(*
_output_shapes	
:
x
PLfc1/bn/gamma/readIdentityPLfc1/bn/gamma*
T0*!
_class
loc:@PLfc1/bn/gamma*
_output_shapes	
:
q
'PLfc1/bn/moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:

PLfc1/bn/moments/meanMeanPLfc1/BiasAdd'PLfc1/bn/moments/mean/reduction_indices*
T0*
_output_shapes
:	*
	keep_dims(*

Tidx0
n
PLfc1/bn/moments/StopGradientStopGradientPLfc1/bn/moments/mean*
T0*
_output_shapes
:	

"PLfc1/bn/moments/SquaredDifferenceSquaredDifferencePLfc1/BiasAddPLfc1/bn/moments/StopGradient*
_output_shapes
:	*
T0
u
+PLfc1/bn/moments/variance/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:
¹
PLfc1/bn/moments/varianceMean"PLfc1/bn/moments/SquaredDifference+PLfc1/bn/moments/variance/reduction_indices*
_output_shapes
:	*
	keep_dims(*

Tidx0*
T0
w
PLfc1/bn/moments/SqueezeSqueezePLfc1/bn/moments/mean*
_output_shapes	
:*
squeeze_dims
 *
T0
}
PLfc1/bn/moments/Squeeze_1SqueezePLfc1/bn/moments/variance*
_output_shapes	
:*
squeeze_dims
 *
T0
_
PLfc1/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

[
PLfc1/bn/cond/switch_tIdentityPLfc1/bn/cond/Switch:1*
_output_shapes
: *
T0

Y
PLfc1/bn/cond/switch_fIdentityPLfc1/bn/cond/Switch*
T0
*
_output_shapes
: 
Q
PLfc1/bn/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
ê
LPLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes	
:
÷
:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
dtype0*
_output_shapes	
:*
shared_name *M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:
ó
APLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverageLPLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
ü
?PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
T0*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
î
NPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes	
:
û
<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes	
:*
shared_name 
û
CPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssign<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverageNPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
_output_shapes	
:*
use_locking(*
T0*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(

APLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentity<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:*
T0*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage

,PLfc1/bn/cond/ExponentialMovingAverage/decayConst^PLfc1/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 

<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^PLfc1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Î
:PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x,PLfc1/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
ù
<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubEPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1GPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_output_shapes	
:
¡
CPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch?PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/readPLfc1/bn/cond/pred_id*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::*
T0
Ú
EPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1SwitchPLfc1/bn/moments/SqueezePLfc1/bn/cond/pred_id*
T0*+
_class!
loc:@PLfc1/bn/moments/Squeeze*"
_output_shapes
::
á
:PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes	
:
È
6PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub?PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1:PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:*
use_locking( *
T0

=PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAveragePLfc1/bn/cond/pred_id*
T0*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^PLfc1/bn/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
Ò
<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x,PLfc1/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
ÿ
>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubGPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1IPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_output_shapes	
:
§
EPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchAPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/readPLfc1/bn/cond/pred_id*"
_output_shapes
::*
T0*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage
à
GPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1SwitchPLfc1/bn/moments/Squeeze_1PLfc1/bn/cond/pred_id*"
_output_shapes
::*
T0*-
_class#
!loc:@PLfc1/bn/moments/Squeeze_1
ç
<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
_output_shapes	
:*
T0
Ð
8PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubAPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
_output_shapes	
:*
use_locking( *
T0*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage

?PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitch<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAveragePLfc1/bn/cond/pred_id*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::*
T0
¢
&PLfc1/bn/cond/ExponentialMovingAverageNoOp7^PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg9^PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
¹
 PLfc1/bn/cond/control_dependencyIdentityPLfc1/bn/cond/switch_t'^PLfc1/bn/cond/ExponentialMovingAverage*)
_class
loc:@PLfc1/bn/cond/switch_t*
_output_shapes
: *
T0

3
PLfc1/bn/cond/NoOpNoOp^PLfc1/bn/cond/switch_f
§
"PLfc1/bn/cond/control_dependency_1IdentityPLfc1/bn/cond/switch_f^PLfc1/bn/cond/NoOp*
T0
*)
_class
loc:@PLfc1/bn/cond/switch_f*
_output_shapes
: 

PLfc1/bn/cond/MergeMerge"PLfc1/bn/cond/control_dependency_1 PLfc1/bn/cond/control_dependency*
N*
_output_shapes
: : *
T0

a
PLfc1/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 
_
PLfc1/bn/cond_1/switch_tIdentityPLfc1/bn/cond_1/Switch:1*
_output_shapes
: *
T0

]
PLfc1/bn/cond_1/switch_fIdentityPLfc1/bn/cond_1/Switch*
_output_shapes
: *
T0

S
PLfc1/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 

PLfc1/bn/cond_1/IdentityIdentity!PLfc1/bn/cond_1/Identity/Switch:1^PLfc1/bn/cond/Merge*
T0*
_output_shapes	
:
¶
PLfc1/bn/cond_1/Identity/SwitchSwitchPLfc1/bn/moments/SqueezePLfc1/bn/cond_1/pred_id*
T0*+
_class!
loc:@PLfc1/bn/moments/Squeeze*"
_output_shapes
::

PLfc1/bn/cond_1/Identity_1Identity#PLfc1/bn/cond_1/Identity_1/Switch:1^PLfc1/bn/cond/Merge*
T0*
_output_shapes	
:
¼
!PLfc1/bn/cond_1/Identity_1/SwitchSwitchPLfc1/bn/moments/Squeeze_1PLfc1/bn/cond_1/pred_id*
T0*-
_class#
!loc:@PLfc1/bn/moments/Squeeze_1*"
_output_shapes
::
ø
PLfc1/bn/cond_1/Switch_1Switch?PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/readPLfc1/bn/cond_1/pred_id*"
_output_shapes
::*
T0*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage
ü
PLfc1/bn/cond_1/Switch_2SwitchAPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/readPLfc1/bn/cond_1/pred_id*
T0*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::

PLfc1/bn/cond_1/MergeMergePLfc1/bn/cond_1/Switch_1PLfc1/bn/cond_1/Identity*
T0*
N*
_output_shapes
	:: 

PLfc1/bn/cond_1/Merge_1MergePLfc1/bn/cond_1/Switch_2PLfc1/bn/cond_1/Identity_1*
T0*
N*
_output_shapes
	:: 
]
PLfc1/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
v
PLfc1/bn/batchnorm/addAddPLfc1/bn/cond_1/Merge_1PLfc1/bn/batchnorm/add/y*
_output_shapes	
:*
T0
_
PLfc1/bn/batchnorm/RsqrtRsqrtPLfc1/bn/batchnorm/add*
T0*
_output_shapes	
:
r
PLfc1/bn/batchnorm/mulMulPLfc1/bn/batchnorm/RsqrtPLfc1/bn/gamma/read*
_output_shapes	
:*
T0
p
PLfc1/bn/batchnorm/mul_1MulPLfc1/BiasAddPLfc1/bn/batchnorm/mul*
T0*
_output_shapes
:	
t
PLfc1/bn/batchnorm/mul_2MulPLfc1/bn/cond_1/MergePLfc1/bn/batchnorm/mul*
T0*
_output_shapes	
:
q
PLfc1/bn/batchnorm/subSubPLfc1/bn/beta/readPLfc1/bn/batchnorm/mul_2*
T0*
_output_shapes	
:
{
PLfc1/bn/batchnorm/add_1AddPLfc1/bn/batchnorm/mul_1PLfc1/bn/batchnorm/sub*
_output_shapes
:	*
T0
V

PLfc1/ReluReluPLfc1/bn/batchnorm/add_1*
_output_shapes
:	*
T0
\
PLdp1/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

U
PLdp1/cond/switch_tIdentityPLdp1/cond/Switch:1*
_output_shapes
: *
T0

S
PLdp1/cond/switch_fIdentityPLdp1/cond/Switch*
T0
*
_output_shapes
: 
N
PLdp1/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

r
PLdp1/cond/dropout/rateConst^PLdp1/cond/switch_t*
dtype0*
_output_shapes
: *
valueB
 *ÍÌÌ>

PLdp1/cond/dropout/ShapeConst^PLdp1/cond/switch_t*
valueB"      *
dtype0*
_output_shapes
:

%PLdp1/cond/dropout/random_uniform/minConst^PLdp1/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

%PLdp1/cond/dropout/random_uniform/maxConst^PLdp1/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
ª
/PLdp1/cond/dropout/random_uniform/RandomUniformRandomUniformPLdp1/cond/dropout/Shape*
T0*
dtype0*
_output_shapes
:	*
seed2 *

seed 

%PLdp1/cond/dropout/random_uniform/subSub%PLdp1/cond/dropout/random_uniform/max%PLdp1/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
®
%PLdp1/cond/dropout/random_uniform/mulMul/PLdp1/cond/dropout/random_uniform/RandomUniform%PLdp1/cond/dropout/random_uniform/sub*
_output_shapes
:	*
T0
 
!PLdp1/cond/dropout/random_uniformAdd%PLdp1/cond/dropout/random_uniform/mul%PLdp1/cond/dropout/random_uniform/min*
T0*
_output_shapes
:	
s
PLdp1/cond/dropout/sub/xConst^PLdp1/cond/switch_t*
_output_shapes
: *
valueB
 *  ?*
dtype0
q
PLdp1/cond/dropout/subSubPLdp1/cond/dropout/sub/xPLdp1/cond/dropout/rate*
_output_shapes
: *
T0
w
PLdp1/cond/dropout/truediv/xConst^PLdp1/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
|
PLdp1/cond/dropout/truedivRealDivPLdp1/cond/dropout/truediv/xPLdp1/cond/dropout/sub*
T0*
_output_shapes
: 

PLdp1/cond/dropout/GreaterEqualGreaterEqual!PLdp1/cond/dropout/random_uniformPLdp1/cond/dropout/rate*
_output_shapes
:	*
T0

PLdp1/cond/dropout/mulMulPLdp1/cond/dropout/mul/Switch:1PLdp1/cond/dropout/truediv*
_output_shapes
:	*
T0

PLdp1/cond/dropout/mul/SwitchSwitch
PLfc1/ReluPLdp1/cond/pred_id*
T0*
_class
loc:@PLfc1/Relu**
_output_shapes
:	:	

PLdp1/cond/dropout/CastCastPLdp1/cond/dropout/GreaterEqual*

SrcT0
*
Truncate( *
_output_shapes
:	*

DstT0
z
PLdp1/cond/dropout/mul_1MulPLdp1/cond/dropout/mulPLdp1/cond/dropout/Cast*
T0*
_output_shapes
:	

PLdp1/cond/Switch_1Switch
PLfc1/ReluPLdp1/cond/pred_id*
_class
loc:@PLfc1/Relu**
_output_shapes
:	:	*
T0
}
PLdp1/cond/MergeMergePLdp1/cond/Switch_1PLdp1/cond/dropout/mul_1*
T0*
N*!
_output_shapes
:	: 
¡
.PLfc2/weights/Initializer/random_uniform/shapeConst*
valueB"      * 
_class
loc:@PLfc2/weights*
dtype0*
_output_shapes
:

,PLfc2/weights/Initializer/random_uniform/minConst*
valueB
 *qÄ¾* 
_class
loc:@PLfc2/weights*
dtype0*
_output_shapes
: 

,PLfc2/weights/Initializer/random_uniform/maxConst*
valueB
 *qÄ>* 
_class
loc:@PLfc2/weights*
dtype0*
_output_shapes
: 
ê
6PLfc2/weights/Initializer/random_uniform/RandomUniformRandomUniform.PLfc2/weights/Initializer/random_uniform/shape*
dtype0* 
_output_shapes
:
*

seed *
T0* 
_class
loc:@PLfc2/weights*
seed2 
Ò
,PLfc2/weights/Initializer/random_uniform/subSub,PLfc2/weights/Initializer/random_uniform/max,PLfc2/weights/Initializer/random_uniform/min*
T0* 
_class
loc:@PLfc2/weights*
_output_shapes
: 
æ
,PLfc2/weights/Initializer/random_uniform/mulMul6PLfc2/weights/Initializer/random_uniform/RandomUniform,PLfc2/weights/Initializer/random_uniform/sub* 
_class
loc:@PLfc2/weights* 
_output_shapes
:
*
T0
Ø
(PLfc2/weights/Initializer/random_uniformAdd,PLfc2/weights/Initializer/random_uniform/mul,PLfc2/weights/Initializer/random_uniform/min* 
_class
loc:@PLfc2/weights* 
_output_shapes
:
*
T0
¶
PLfc2/weights
VariableV2"/device:CPU:0*
shared_name * 
_class
loc:@PLfc2/weights*
	container *
shape:
*
dtype0* 
_output_shapes
:

Ü
PLfc2/weights/AssignAssignPLfc2/weights(PLfc2/weights/Initializer/random_uniform"/device:CPU:0*
validate_shape(* 
_output_shapes
:
*
use_locking(*
T0* 
_class
loc:@PLfc2/weights

PLfc2/weights/readIdentityPLfc2/weights"/device:CPU:0*
T0* 
_class
loc:@PLfc2/weights* 
_output_shapes
:

K
PLfc2/L2LossL2LossPLfc2/weights/read*
T0*
_output_shapes
: 
X
PLfc2/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
\
PLfc2/weight_lossMulPLfc2/L2LossPLfc2/weight_loss/y*
T0*
_output_shapes
: 

PLfc2/MatMulMatMulPLdp1/cond/MergePLfc2/weights/read*
T0*
_output_shapes
:	*
transpose_a( *
transpose_b( 

PLfc2/biases/Initializer/ConstConst*
_output_shapes	
:*
valueB*    *
_class
loc:@PLfc2/biases*
dtype0
ª
PLfc2/biases
VariableV2"/device:CPU:0*
dtype0*
_output_shapes	
:*
shared_name *
_class
loc:@PLfc2/biases*
	container *
shape:
Ê
PLfc2/biases/AssignAssignPLfc2/biasesPLfc2/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*
_class
loc:@PLfc2/biases*
validate_shape(*
_output_shapes	
:

PLfc2/biases/readIdentityPLfc2/biases"/device:CPU:0*
T0*
_class
loc:@PLfc2/biases*
_output_shapes	
:
z
PLfc2/BiasAddBiasAddPLfc2/MatMulPLfc2/biases/read*
data_formatNHWC*
_output_shapes
:	*
T0
]
PLfc2/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
{
PLfc2/bn/beta
VariableV2*
shared_name *
dtype0*
_output_shapes	
:*
	container *
shape:
®
PLfc2/bn/beta/AssignAssignPLfc2/bn/betaPLfc2/bn/Const*
_output_shapes	
:*
use_locking(*
T0* 
_class
loc:@PLfc2/bn/beta*
validate_shape(
u
PLfc2/bn/beta/readIdentityPLfc2/bn/beta*
T0* 
_class
loc:@PLfc2/bn/beta*
_output_shapes	
:
_
PLfc2/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:
|
PLfc2/bn/gamma
VariableV2*
dtype0*
_output_shapes	
:*
	container *
shape:*
shared_name 
³
PLfc2/bn/gamma/AssignAssignPLfc2/bn/gammaPLfc2/bn/Const_1*
_output_shapes	
:*
use_locking(*
T0*!
_class
loc:@PLfc2/bn/gamma*
validate_shape(
x
PLfc2/bn/gamma/readIdentityPLfc2/bn/gamma*!
_class
loc:@PLfc2/bn/gamma*
_output_shapes	
:*
T0
q
'PLfc2/bn/moments/mean/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:

PLfc2/bn/moments/meanMeanPLfc2/BiasAdd'PLfc2/bn/moments/mean/reduction_indices*
T0*
_output_shapes
:	*
	keep_dims(*

Tidx0
n
PLfc2/bn/moments/StopGradientStopGradientPLfc2/bn/moments/mean*
T0*
_output_shapes
:	

"PLfc2/bn/moments/SquaredDifferenceSquaredDifferencePLfc2/BiasAddPLfc2/bn/moments/StopGradient*
T0*
_output_shapes
:	
u
+PLfc2/bn/moments/variance/reduction_indicesConst*
valueB: *
dtype0*
_output_shapes
:
¹
PLfc2/bn/moments/varianceMean"PLfc2/bn/moments/SquaredDifference+PLfc2/bn/moments/variance/reduction_indices*
_output_shapes
:	*
	keep_dims(*

Tidx0*
T0
w
PLfc2/bn/moments/SqueezeSqueezePLfc2/bn/moments/mean*
_output_shapes	
:*
squeeze_dims
 *
T0
}
PLfc2/bn/moments/Squeeze_1SqueezePLfc2/bn/moments/variance*
_output_shapes	
:*
squeeze_dims
 *
T0
_
PLfc2/bn/cond/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

[
PLfc2/bn/cond/switch_tIdentityPLfc2/bn/cond/Switch:1*
T0
*
_output_shapes
: 
Y
PLfc2/bn/cond/switch_fIdentityPLfc2/bn/cond/Switch*
_output_shapes
: *
T0

Q
PLfc2/bn/cond/pred_idIdentityPlaceholder_2*
_output_shapes
: *
T0

ê
LPLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosConst*
valueB*    *M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
dtype0*
_output_shapes	
:
÷
:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
_output_shapes	
:*
shared_name *M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
	container *
shape:*
dtype0
ó
APLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverageLPLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
_output_shapes	
:*
use_locking(*
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(
ü
?PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
î
NPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosConst*
_output_shapes	
:*
valueB*    *O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0
û
<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*
	container *
shape:*
dtype0*
_output_shapes	
:*
shared_name 
û
CPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssign<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverageNPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
use_locking(*
T0*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

APLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentity<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:*
T0*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage

,PLfc2/bn/cond/ExponentialMovingAverage/decayConst^PLfc2/bn/cond/switch_t*
_output_shapes
: *
valueB
 *fff?*
dtype0

<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^PLfc2/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Î
:PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x,PLfc2/bn/cond/ExponentialMovingAverage/decay*
T0*
_output_shapes
: 
ù
<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubEPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1GPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
_output_shapes	
:*
T0
¡
CPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch?PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/readPLfc2/bn/cond/pred_id*
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
Ú
EPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1SwitchPLfc2/bn/moments/SqueezePLfc2/bn/cond/pred_id*+
_class!
loc:@PLfc2/bn/moments/Squeeze*"
_output_shapes
::*
T0
á
:PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_output_shapes	
:
È
6PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub?PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1:PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
use_locking( *
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

=PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAveragePLfc2/bn/cond/pred_id*
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^PLfc2/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
Ò
<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x,PLfc2/bn/cond/ExponentialMovingAverage/decay*
_output_shapes
: *
T0
ÿ
>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubGPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1IPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
_output_shapes	
:*
T0
§
EPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchAPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/readPLfc2/bn/cond/pred_id*
T0*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
à
GPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1SwitchPLfc2/bn/moments/Squeeze_1PLfc2/bn/cond/pred_id*"
_output_shapes
::*
T0*-
_class#
!loc:@PLfc2/bn/moments/Squeeze_1
ç
<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_output_shapes	
:
Ð
8PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubAPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
use_locking( *
T0*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

?PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitch<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAveragePLfc2/bn/cond/pred_id*
T0*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
¢
&PLfc2/bn/cond/ExponentialMovingAverageNoOp7^PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg9^PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1
¹
 PLfc2/bn/cond/control_dependencyIdentityPLfc2/bn/cond/switch_t'^PLfc2/bn/cond/ExponentialMovingAverage*
_output_shapes
: *
T0
*)
_class
loc:@PLfc2/bn/cond/switch_t
3
PLfc2/bn/cond/NoOpNoOp^PLfc2/bn/cond/switch_f
§
"PLfc2/bn/cond/control_dependency_1IdentityPLfc2/bn/cond/switch_f^PLfc2/bn/cond/NoOp*
_output_shapes
: *
T0
*)
_class
loc:@PLfc2/bn/cond/switch_f

PLfc2/bn/cond/MergeMerge"PLfc2/bn/cond/control_dependency_1 PLfc2/bn/cond/control_dependency*
N*
_output_shapes
: : *
T0

a
PLfc2/bn/cond_1/SwitchSwitchPlaceholder_2Placeholder_2*
_output_shapes
: : *
T0

_
PLfc2/bn/cond_1/switch_tIdentityPLfc2/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
]
PLfc2/bn/cond_1/switch_fIdentityPLfc2/bn/cond_1/Switch*
_output_shapes
: *
T0

S
PLfc2/bn/cond_1/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 

PLfc2/bn/cond_1/IdentityIdentity!PLfc2/bn/cond_1/Identity/Switch:1^PLfc2/bn/cond/Merge*
T0*
_output_shapes	
:
¶
PLfc2/bn/cond_1/Identity/SwitchSwitchPLfc2/bn/moments/SqueezePLfc2/bn/cond_1/pred_id*"
_output_shapes
::*
T0*+
_class!
loc:@PLfc2/bn/moments/Squeeze

PLfc2/bn/cond_1/Identity_1Identity#PLfc2/bn/cond_1/Identity_1/Switch:1^PLfc2/bn/cond/Merge*
T0*
_output_shapes	
:
¼
!PLfc2/bn/cond_1/Identity_1/SwitchSwitchPLfc2/bn/moments/Squeeze_1PLfc2/bn/cond_1/pred_id*-
_class#
!loc:@PLfc2/bn/moments/Squeeze_1*"
_output_shapes
::*
T0
ø
PLfc2/bn/cond_1/Switch_1Switch?PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/readPLfc2/bn/cond_1/pred_id*
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
ü
PLfc2/bn/cond_1/Switch_2SwitchAPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/readPLfc2/bn/cond_1/pred_id*"
_output_shapes
::*
T0*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage

PLfc2/bn/cond_1/MergeMergePLfc2/bn/cond_1/Switch_1PLfc2/bn/cond_1/Identity*
N*
_output_shapes
	:: *
T0

PLfc2/bn/cond_1/Merge_1MergePLfc2/bn/cond_1/Switch_2PLfc2/bn/cond_1/Identity_1*
N*
_output_shapes
	:: *
T0
]
PLfc2/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
v
PLfc2/bn/batchnorm/addAddPLfc2/bn/cond_1/Merge_1PLfc2/bn/batchnorm/add/y*
_output_shapes	
:*
T0
_
PLfc2/bn/batchnorm/RsqrtRsqrtPLfc2/bn/batchnorm/add*
T0*
_output_shapes	
:
r
PLfc2/bn/batchnorm/mulMulPLfc2/bn/batchnorm/RsqrtPLfc2/bn/gamma/read*
_output_shapes	
:*
T0
p
PLfc2/bn/batchnorm/mul_1MulPLfc2/BiasAddPLfc2/bn/batchnorm/mul*
_output_shapes
:	*
T0
t
PLfc2/bn/batchnorm/mul_2MulPLfc2/bn/cond_1/MergePLfc2/bn/batchnorm/mul*
T0*
_output_shapes	
:
q
PLfc2/bn/batchnorm/subSubPLfc2/bn/beta/readPLfc2/bn/batchnorm/mul_2*
T0*
_output_shapes	
:
{
PLfc2/bn/batchnorm/add_1AddPLfc2/bn/batchnorm/mul_1PLfc2/bn/batchnorm/sub*
_output_shapes
:	*
T0
V

PLfc2/ReluReluPLfc2/bn/batchnorm/add_1*
T0*
_output_shapes
:	
\
PLdp2/cond/SwitchSwitchPlaceholder_2Placeholder_2*
T0
*
_output_shapes
: : 
U
PLdp2/cond/switch_tIdentityPLdp2/cond/Switch:1*
_output_shapes
: *
T0

S
PLdp2/cond/switch_fIdentityPLdp2/cond/Switch*
_output_shapes
: *
T0

N
PLdp2/cond/pred_idIdentityPlaceholder_2*
T0
*
_output_shapes
: 
r
PLdp2/cond/dropout/rateConst^PLdp2/cond/switch_t*
_output_shapes
: *
valueB
 *ÍÌÌ>*
dtype0

PLdp2/cond/dropout/ShapeConst^PLdp2/cond/switch_t*
valueB"      *
dtype0*
_output_shapes
:

%PLdp2/cond/dropout/random_uniform/minConst^PLdp2/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

%PLdp2/cond/dropout/random_uniform/maxConst^PLdp2/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
ª
/PLdp2/cond/dropout/random_uniform/RandomUniformRandomUniformPLdp2/cond/dropout/Shape*
dtype0*
_output_shapes
:	*
seed2 *

seed *
T0

%PLdp2/cond/dropout/random_uniform/subSub%PLdp2/cond/dropout/random_uniform/max%PLdp2/cond/dropout/random_uniform/min*
_output_shapes
: *
T0
®
%PLdp2/cond/dropout/random_uniform/mulMul/PLdp2/cond/dropout/random_uniform/RandomUniform%PLdp2/cond/dropout/random_uniform/sub*
T0*
_output_shapes
:	
 
!PLdp2/cond/dropout/random_uniformAdd%PLdp2/cond/dropout/random_uniform/mul%PLdp2/cond/dropout/random_uniform/min*
T0*
_output_shapes
:	
s
PLdp2/cond/dropout/sub/xConst^PLdp2/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
q
PLdp2/cond/dropout/subSubPLdp2/cond/dropout/sub/xPLdp2/cond/dropout/rate*
T0*
_output_shapes
: 
w
PLdp2/cond/dropout/truediv/xConst^PLdp2/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
|
PLdp2/cond/dropout/truedivRealDivPLdp2/cond/dropout/truediv/xPLdp2/cond/dropout/sub*
T0*
_output_shapes
: 

PLdp2/cond/dropout/GreaterEqualGreaterEqual!PLdp2/cond/dropout/random_uniformPLdp2/cond/dropout/rate*
_output_shapes
:	*
T0

PLdp2/cond/dropout/mulMulPLdp2/cond/dropout/mul/Switch:1PLdp2/cond/dropout/truediv*
T0*
_output_shapes
:	

PLdp2/cond/dropout/mul/SwitchSwitch
PLfc2/ReluPLdp2/cond/pred_id*
_class
loc:@PLfc2/Relu**
_output_shapes
:	:	*
T0

PLdp2/cond/dropout/CastCastPLdp2/cond/dropout/GreaterEqual*

SrcT0
*
Truncate( *
_output_shapes
:	*

DstT0
z
PLdp2/cond/dropout/mul_1MulPLdp2/cond/dropout/mulPLdp2/cond/dropout/Cast*
_output_shapes
:	*
T0

PLdp2/cond/Switch_1Switch
PLfc2/ReluPLdp2/cond/pred_id*
T0*
_class
loc:@PLfc2/Relu**
_output_shapes
:	:	
}
PLdp2/cond/MergeMergePLdp2/cond/Switch_1PLdp2/cond/dropout/mul_1*
T0*
N*!
_output_shapes
:	: 
¡
.fc3PL/weights/Initializer/random_uniform/shapeConst*
_output_shapes
:*
valueB"      * 
_class
loc:@fc3PL/weights*
dtype0

,fc3PL/weights/Initializer/random_uniform/minConst*
_output_shapes
: *
valueB
 *ý[¾* 
_class
loc:@fc3PL/weights*
dtype0

,fc3PL/weights/Initializer/random_uniform/maxConst*
valueB
 *ý[>* 
_class
loc:@fc3PL/weights*
dtype0*
_output_shapes
: 
é
6fc3PL/weights/Initializer/random_uniform/RandomUniformRandomUniform.fc3PL/weights/Initializer/random_uniform/shape*
dtype0*
_output_shapes
:	*

seed *
T0* 
_class
loc:@fc3PL/weights*
seed2 
Ò
,fc3PL/weights/Initializer/random_uniform/subSub,fc3PL/weights/Initializer/random_uniform/max,fc3PL/weights/Initializer/random_uniform/min*
T0* 
_class
loc:@fc3PL/weights*
_output_shapes
: 
å
,fc3PL/weights/Initializer/random_uniform/mulMul6fc3PL/weights/Initializer/random_uniform/RandomUniform,fc3PL/weights/Initializer/random_uniform/sub*
_output_shapes
:	*
T0* 
_class
loc:@fc3PL/weights
×
(fc3PL/weights/Initializer/random_uniformAdd,fc3PL/weights/Initializer/random_uniform/mul,fc3PL/weights/Initializer/random_uniform/min*
T0* 
_class
loc:@fc3PL/weights*
_output_shapes
:	
´
fc3PL/weights
VariableV2"/device:CPU:0*
shared_name * 
_class
loc:@fc3PL/weights*
	container *
shape:	*
dtype0*
_output_shapes
:	
Û
fc3PL/weights/AssignAssignfc3PL/weights(fc3PL/weights/Initializer/random_uniform"/device:CPU:0*
_output_shapes
:	*
use_locking(*
T0* 
_class
loc:@fc3PL/weights*
validate_shape(

fc3PL/weights/readIdentityfc3PL/weights"/device:CPU:0*
T0* 
_class
loc:@fc3PL/weights*
_output_shapes
:	
K
fc3PL/L2LossL2Lossfc3PL/weights/read*
_output_shapes
: *
T0
X
fc3PL/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
\
fc3PL/weight_lossMulfc3PL/L2Lossfc3PL/weight_loss/y*
_output_shapes
: *
T0

fc3PL/MatMulMatMulPLdp2/cond/Mergefc3PL/weights/read*
T0*
_output_shapes

:*
transpose_a( *
transpose_b( 

fc3PL/biases/Initializer/ConstConst*
valueB*    *
_class
loc:@fc3PL/biases*
dtype0*
_output_shapes
:
¨
fc3PL/biases
VariableV2"/device:CPU:0*
_output_shapes
:*
shared_name *
_class
loc:@fc3PL/biases*
	container *
shape:*
dtype0
É
fc3PL/biases/AssignAssignfc3PL/biasesfc3PL/biases/Initializer/Const"/device:CPU:0*
use_locking(*
T0*
_class
loc:@fc3PL/biases*
validate_shape(*
_output_shapes
:

fc3PL/biases/readIdentityfc3PL/biases"/device:CPU:0*
T0*
_class
loc:@fc3PL/biases*
_output_shapes
:
y
fc3PL/BiasAddBiasAddfc3PL/MatMulfc3PL/biases/read*
T0*
data_formatNHWC*
_output_shapes

:
L
	Softmax_2Softmaxfc3PL/BiasAdd*
T0*
_output_shapes

:
Y
save/filename/inputConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
n
save/filenamePlaceholderWithDefaultsave/filename/input*
dtype0*
_output_shapes
: *
shape: 
e

save/ConstPlaceholderWithDefaultsave/filename*
dtype0*
_output_shapes
: *
shape: 
Ã$
save/SaveV2/tensor_namesConst*
_output_shapes
:\*ö#
valueì#Bé#\BPLfc1/biasesB:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc1/bn/betaBPLfc1/bn/gammaBPLfc1/weightsBPLfc2/biasesB:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc2/bn/betaBPLfc2/bn/gammaBPLfc2/weightsBaggPL/biasesB:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverageB<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageBaggPL/bn/betaBaggPL/bn/gammaBaggPL/weightsBfc3PL/biasesBfc3PL/weightsBgapnet00/biasesBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet00/weightsBgapnet01PL/biasesBgapnet01PL/bn/betaBgapnet01PL/bn/gammaBDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01PL/weightsBgapnet10/biasesBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/weightsBgapnet11PL/biasesBgapnet11PL/bn/betaBgapnet11PL/bn/gammaBDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11PL/weightsBlayerfilter0PL/BiasAdd/biasesBlayerfilter0PL_edgefea_0/biasesB layerfilter0PL_edgefea_0/bn/betaB!layerfilter0PL_edgefea_0/bn/gammaB`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB layerfilter0PL_edgefea_0/weightsB*layerfilter0PL_neib_att_conv_head_0/biasesB+layerfilter0PL_neib_att_conv_head_0/bn/betaB,layerfilter0PL_neib_att_conv_head_0/bn/gammaBvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_neib_att_conv_head_0/weightsB)layerfilter0PL_newfea_conv_head_0/bn/betaB*layerfilter0PL_newfea_conv_head_0/bn/gammaBrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0PL_newfea_conv_head_0/weightsB*layerfilter0PL_self_att_conv_head_0/biasesB+layerfilter0PL_self_att_conv_head_0/bn/betaB,layerfilter0PL_self_att_conv_head_0/bn/gammaBvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_self_att_conv_head_0/weightsBlayerfilter1PL/BiasAdd/biasesBlayerfilter1PL_edgefea_0/biasesB layerfilter1PL_edgefea_0/bn/betaB!layerfilter1PL_edgefea_0/bn/gammaB`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB layerfilter1PL_edgefea_0/weightsB*layerfilter1PL_neib_att_conv_head_0/biasesB+layerfilter1PL_neib_att_conv_head_0/bn/betaB,layerfilter1PL_neib_att_conv_head_0/bn/gammaBvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_neib_att_conv_head_0/weightsB)layerfilter1PL_newfea_conv_head_0/bn/betaB*layerfilter1PL_newfea_conv_head_0/bn/gammaBrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1PL_newfea_conv_head_0/weightsB*layerfilter1PL_self_att_conv_head_0/biasesB+layerfilter1PL_self_att_conv_head_0/bn/betaB,layerfilter1PL_self_att_conv_head_0/bn/gammaBvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_self_att_conv_head_0/weights*
dtype0

save/SaveV2/shape_and_slicesConst*Í
valueÃBÀ\B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:\
¦%
save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicesPLfc1/biases:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAveragePLfc1/bn/betaPLfc1/bn/gammaPLfc1/weightsPLfc2/biases:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAveragePLfc2/bn/betaPLfc2/bn/gammaPLfc2/weightsaggPL/biases:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageaggPL/bn/betaaggPL/bn/gammaaggPL/weightsfc3PL/biasesfc3PL/weightsgapnet00/biasesgapnet00/bn/betagapnet00/bn/gamma@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet00/weightsgapnet01PL/biasesgapnet01PL/bn/betagapnet01PL/bn/gammaDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverageFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet01PL/weightsgapnet10/biasesgapnet10/bn/betagapnet10/bn/gamma@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet10/weightsgapnet11PL/biasesgapnet11PL/bn/betagapnet11PL/bn/gammaDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverageFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet11PL/weightslayerfilter0PL/BiasAdd/biaseslayerfilter0PL_edgefea_0/biases layerfilter0PL_edgefea_0/bn/beta!layerfilter0PL_edgefea_0/bn/gamma`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage layerfilter0PL_edgefea_0/weights*layerfilter0PL_neib_att_conv_head_0/biases+layerfilter0PL_neib_att_conv_head_0/bn/beta,layerfilter0PL_neib_att_conv_head_0/bn/gammavlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter0PL_neib_att_conv_head_0/weights)layerfilter0PL_newfea_conv_head_0/bn/beta*layerfilter0PL_newfea_conv_head_0/bn/gammarlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter0PL_newfea_conv_head_0/weights*layerfilter0PL_self_att_conv_head_0/biases+layerfilter0PL_self_att_conv_head_0/bn/beta,layerfilter0PL_self_att_conv_head_0/bn/gammavlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter0PL_self_att_conv_head_0/weightslayerfilter1PL/BiasAdd/biaseslayerfilter1PL_edgefea_0/biases layerfilter1PL_edgefea_0/bn/beta!layerfilter1PL_edgefea_0/bn/gamma`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage layerfilter1PL_edgefea_0/weights*layerfilter1PL_neib_att_conv_head_0/biases+layerfilter1PL_neib_att_conv_head_0/bn/beta,layerfilter1PL_neib_att_conv_head_0/bn/gammavlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter1PL_neib_att_conv_head_0/weights)layerfilter1PL_newfea_conv_head_0/bn/beta*layerfilter1PL_newfea_conv_head_0/bn/gammarlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1PL_newfea_conv_head_0/weights*layerfilter1PL_self_att_conv_head_0/biases+layerfilter1PL_self_att_conv_head_0/bn/beta,layerfilter1PL_self_att_conv_head_0/bn/gammavlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter1PL_self_att_conv_head_0/weights*j
dtypes`
^2\
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 
Õ$
save/RestoreV2/tensor_namesConst"/device:CPU:0*ö#
valueì#Bé#\BPLfc1/biasesB:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc1/bn/betaBPLfc1/bn/gammaBPLfc1/weightsBPLfc2/biasesB:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc2/bn/betaBPLfc2/bn/gammaBPLfc2/weightsBaggPL/biasesB:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverageB<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageBaggPL/bn/betaBaggPL/bn/gammaBaggPL/weightsBfc3PL/biasesBfc3PL/weightsBgapnet00/biasesBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet00/weightsBgapnet01PL/biasesBgapnet01PL/bn/betaBgapnet01PL/bn/gammaBDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01PL/weightsBgapnet10/biasesBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/weightsBgapnet11PL/biasesBgapnet11PL/bn/betaBgapnet11PL/bn/gammaBDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11PL/weightsBlayerfilter0PL/BiasAdd/biasesBlayerfilter0PL_edgefea_0/biasesB layerfilter0PL_edgefea_0/bn/betaB!layerfilter0PL_edgefea_0/bn/gammaB`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB layerfilter0PL_edgefea_0/weightsB*layerfilter0PL_neib_att_conv_head_0/biasesB+layerfilter0PL_neib_att_conv_head_0/bn/betaB,layerfilter0PL_neib_att_conv_head_0/bn/gammaBvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_neib_att_conv_head_0/weightsB)layerfilter0PL_newfea_conv_head_0/bn/betaB*layerfilter0PL_newfea_conv_head_0/bn/gammaBrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0PL_newfea_conv_head_0/weightsB*layerfilter0PL_self_att_conv_head_0/biasesB+layerfilter0PL_self_att_conv_head_0/bn/betaB,layerfilter0PL_self_att_conv_head_0/bn/gammaBvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_self_att_conv_head_0/weightsBlayerfilter1PL/BiasAdd/biasesBlayerfilter1PL_edgefea_0/biasesB layerfilter1PL_edgefea_0/bn/betaB!layerfilter1PL_edgefea_0/bn/gammaB`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB layerfilter1PL_edgefea_0/weightsB*layerfilter1PL_neib_att_conv_head_0/biasesB+layerfilter1PL_neib_att_conv_head_0/bn/betaB,layerfilter1PL_neib_att_conv_head_0/bn/gammaBvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_neib_att_conv_head_0/weightsB)layerfilter1PL_newfea_conv_head_0/bn/betaB*layerfilter1PL_newfea_conv_head_0/bn/gammaBrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1PL_newfea_conv_head_0/weightsB*layerfilter1PL_self_att_conv_head_0/biasesB+layerfilter1PL_self_att_conv_head_0/bn/betaB,layerfilter1PL_self_att_conv_head_0/bn/gammaBvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
:\
°
save/RestoreV2/shape_and_slicesConst"/device:CPU:0*Í
valueÃBÀ\B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:\
é
save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*
_output_shapesó
ð::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*j
dtypes`
^2\
²
save/AssignAssignPLfc1/biasessave/RestoreV2"/device:CPU:0*
use_locking(*
T0*
_class
loc:@PLfc1/biases*
validate_shape(*
_output_shapes	
:

save/Assign_1Assign:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:1*
use_locking(*
T0*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save/Assign_2Assign<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:2*
use_locking(*
T0*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
©
save/Assign_3AssignPLfc1/bn/betasave/RestoreV2:3*
use_locking(*
T0* 
_class
loc:@PLfc1/bn/beta*
validate_shape(*
_output_shapes	
:
«
save/Assign_4AssignPLfc1/bn/gammasave/RestoreV2:4*
_output_shapes	
:*
use_locking(*
T0*!
_class
loc:@PLfc1/bn/gamma*
validate_shape(
½
save/Assign_5AssignPLfc1/weightssave/RestoreV2:5"/device:CPU:0*
use_locking(*
T0* 
_class
loc:@PLfc1/weights*
validate_shape(* 
_output_shapes
:

¶
save/Assign_6AssignPLfc2/biasessave/RestoreV2:6"/device:CPU:0*
_output_shapes	
:*
use_locking(*
T0*
_class
loc:@PLfc2/biases*
validate_shape(

save/Assign_7Assign:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:7*
use_locking(*
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save/Assign_8Assign<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:8*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
©
save/Assign_9AssignPLfc2/bn/betasave/RestoreV2:9*
_output_shapes	
:*
use_locking(*
T0* 
_class
loc:@PLfc2/bn/beta*
validate_shape(
­
save/Assign_10AssignPLfc2/bn/gammasave/RestoreV2:10*!
_class
loc:@PLfc2/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
¿
save/Assign_11AssignPLfc2/weightssave/RestoreV2:11"/device:CPU:0*
use_locking(*
T0* 
_class
loc:@PLfc2/weights*
validate_shape(* 
_output_shapes
:

¸
save/Assign_12AssignaggPL/biasessave/RestoreV2:12"/device:CPU:0*
use_locking(*
T0*
_class
loc:@aggPL/biases*
validate_shape(*
_output_shapes	
:

save/Assign_13Assign:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:13*
use_locking(*
T0*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save/Assign_14Assign<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:14*
use_locking(*
T0*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
«
save/Assign_15AssignaggPL/bn/betasave/RestoreV2:15* 
_class
loc:@aggPL/bn/beta*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
­
save/Assign_16AssignaggPL/bn/gammasave/RestoreV2:16*!
_class
loc:@aggPL/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
Ç
save/Assign_17AssignaggPL/weightssave/RestoreV2:17"/device:CPU:0*
use_locking(*
T0* 
_class
loc:@aggPL/weights*
validate_shape(*(
_output_shapes
:
·
save/Assign_18Assignfc3PL/biasessave/RestoreV2:18"/device:CPU:0*
use_locking(*
T0*
_class
loc:@fc3PL/biases*
validate_shape(*
_output_shapes
:
¾
save/Assign_19Assignfc3PL/weightssave/RestoreV2:19"/device:CPU:0*
use_locking(*
T0* 
_class
loc:@fc3PL/weights*
validate_shape(*
_output_shapes
:	
½
save/Assign_20Assigngapnet00/biasessave/RestoreV2:20"/device:CPU:0*
use_locking(*
T0*"
_class
loc:@gapnet00/biases*
validate_shape(*
_output_shapes
:@
°
save/Assign_21Assigngapnet00/bn/betasave/RestoreV2:21*#
_class
loc:@gapnet00/bn/beta*
validate_shape(*
_output_shapes
:@*
use_locking(*
T0
²
save/Assign_22Assigngapnet00/bn/gammasave/RestoreV2:22*
_output_shapes
:@*
use_locking(*
T0*$
_class
loc:@gapnet00/bn/gamma*
validate_shape(

save/Assign_23Assign@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:23*
use_locking(*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@

save/Assign_24AssignBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:24*
_output_shapes
:@*
use_locking(*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
Ë
save/Assign_25Assigngapnet00/weightssave/RestoreV2:25"/device:CPU:0*#
_class
loc:@gapnet00/weights*
validate_shape(*&
_output_shapes
:@*
use_locking(*
T0
Â
save/Assign_26Assigngapnet01PL/biasessave/RestoreV2:26"/device:CPU:0*$
_class
loc:@gapnet01PL/biases*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
µ
save/Assign_27Assigngapnet01PL/bn/betasave/RestoreV2:27*%
_class
loc:@gapnet01PL/bn/beta*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
·
save/Assign_28Assigngapnet01PL/bn/gammasave/RestoreV2:28*
use_locking(*
T0*&
_class
loc:@gapnet01PL/bn/gamma*
validate_shape(*
_output_shapes	
:

save/Assign_29AssignDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:29*
_output_shapes	
:*
use_locking(*
T0*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(

save/Assign_30AssignFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:30*
_output_shapes	
:*
use_locking(*
T0*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
Ð
save/Assign_31Assigngapnet01PL/weightssave/RestoreV2:31"/device:CPU:0*'
_output_shapes
:@*
use_locking(*
T0*%
_class
loc:@gapnet01PL/weights*
validate_shape(
¾
save/Assign_32Assigngapnet10/biasessave/RestoreV2:32"/device:CPU:0*"
_class
loc:@gapnet10/biases*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
±
save/Assign_33Assigngapnet10/bn/betasave/RestoreV2:33*
use_locking(*
T0*#
_class
loc:@gapnet10/bn/beta*
validate_shape(*
_output_shapes	
:
³
save/Assign_34Assigngapnet10/bn/gammasave/RestoreV2:34*
use_locking(*
T0*$
_class
loc:@gapnet10/bn/gamma*
validate_shape(*
_output_shapes	
:

save/Assign_35Assign@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:35*
_output_shapes	
:*
use_locking(*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(

save/Assign_36AssignBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:36*
use_locking(*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
Ì
save/Assign_37Assigngapnet10/weightssave/RestoreV2:37"/device:CPU:0*
use_locking(*
T0*#
_class
loc:@gapnet10/weights*
validate_shape(*'
_output_shapes
:M
Â
save/Assign_38Assigngapnet11PL/biasessave/RestoreV2:38"/device:CPU:0*
use_locking(*
T0*$
_class
loc:@gapnet11PL/biases*
validate_shape(*
_output_shapes	
:
µ
save/Assign_39Assigngapnet11PL/bn/betasave/RestoreV2:39*
use_locking(*
T0*%
_class
loc:@gapnet11PL/bn/beta*
validate_shape(*
_output_shapes	
:
·
save/Assign_40Assigngapnet11PL/bn/gammasave/RestoreV2:40*&
_class
loc:@gapnet11PL/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

save/Assign_41AssignDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:41*
use_locking(*
T0*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save/Assign_42AssignFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:42*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
Ñ
save/Assign_43Assigngapnet11PL/weightssave/RestoreV2:43"/device:CPU:0*%
_class
loc:@gapnet11PL/weights*
validate_shape(*(
_output_shapes
:*
use_locking(*
T0
Ê
save/Assign_44Assignlayerfilter0PL/BiasAdd/biasessave/RestoreV2:44*
_output_shapes
:*
use_locking(*
T0*0
_class&
$"loc:@layerfilter0PL/BiasAdd/biases*
validate_shape(
Ý
save/Assign_45Assignlayerfilter0PL_edgefea_0/biasessave/RestoreV2:45"/device:CPU:0*
use_locking(*
T0*2
_class(
&$loc:@layerfilter0PL_edgefea_0/biases*
validate_shape(*
_output_shapes
:
Ð
save/Assign_46Assign layerfilter0PL_edgefea_0/bn/betasave/RestoreV2:46*
use_locking(*
T0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/bn/beta*
validate_shape(*
_output_shapes
:
Ò
save/Assign_47Assign!layerfilter0PL_edgefea_0/bn/gammasave/RestoreV2:47*
use_locking(*
T0*4
_class*
(&loc:@layerfilter0PL_edgefea_0/bn/gamma*
validate_shape(*
_output_shapes
:
Ð
save/Assign_48Assign`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:48*
_output_shapes
:*
use_locking(*
T0*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(
Ô
save/Assign_49Assignblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:49*
use_locking(*
T0*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ë
save/Assign_50Assign layerfilter0PL_edgefea_0/weightssave/RestoreV2:50"/device:CPU:0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
validate_shape(*&
_output_shapes
:*
use_locking(*
T0
ó
save/Assign_51Assign*layerfilter0PL_neib_att_conv_head_0/biasessave/RestoreV2:51"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_neib_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:
æ
save/Assign_52Assign+layerfilter0PL_neib_att_conv_head_0/bn/betasave/RestoreV2:52*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
è
save/Assign_53Assign,layerfilter0PL_neib_att_conv_head_0/bn/gammasave/RestoreV2:53*
use_locking(*
T0*?
_class5
31loc:@layerfilter0PL_neib_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
ý
save/Assign_54Assignvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:54*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:*
use_locking(*
T0

save/Assign_55Assignxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:55*
use_locking(*
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:

save/Assign_56Assign+layerfilter0PL_neib_att_conv_head_0/weightssave/RestoreV2:56"/device:CPU:0*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
validate_shape(*&
_output_shapes
:
â
save/Assign_57Assign)layerfilter0PL_newfea_conv_head_0/bn/betasave/RestoreV2:57*
use_locking(*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
ä
save/Assign_58Assign*layerfilter0PL_newfea_conv_head_0/bn/gammasave/RestoreV2:58*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_newfea_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
õ
save/Assign_59Assignrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:59*
use_locking(*
T0*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ù
save/Assign_60Assigntlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:60*
use_locking(*
T0*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ý
save/Assign_61Assign)layerfilter0PL_newfea_conv_head_0/weightssave/RestoreV2:61"/device:CPU:0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
validate_shape(*&
_output_shapes
:*
use_locking(*
T0
ó
save/Assign_62Assign*layerfilter0PL_self_att_conv_head_0/biasessave/RestoreV2:62"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_self_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:
æ
save/Assign_63Assign+layerfilter0PL_self_att_conv_head_0/bn/betasave/RestoreV2:63*
_output_shapes
:*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/bn/beta*
validate_shape(
è
save/Assign_64Assign,layerfilter0PL_self_att_conv_head_0/bn/gammasave/RestoreV2:64*
use_locking(*
T0*?
_class5
31loc:@layerfilter0PL_self_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
ý
save/Assign_65Assignvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:65*
_output_shapes
:*
use_locking(*
T0*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(

save/Assign_66Assignxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:66*
use_locking(*
T0*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:

save/Assign_67Assign+layerfilter0PL_self_att_conv_head_0/weightssave/RestoreV2:67"/device:CPU:0*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
validate_shape(*&
_output_shapes
:
Ê
save/Assign_68Assignlayerfilter1PL/BiasAdd/biasessave/RestoreV2:68*
use_locking(*
T0*0
_class&
$"loc:@layerfilter1PL/BiasAdd/biases*
validate_shape(*
_output_shapes
:@
Ý
save/Assign_69Assignlayerfilter1PL_edgefea_0/biasessave/RestoreV2:69"/device:CPU:0*
use_locking(*
T0*2
_class(
&$loc:@layerfilter1PL_edgefea_0/biases*
validate_shape(*
_output_shapes
:@
Ð
save/Assign_70Assign layerfilter1PL_edgefea_0/bn/betasave/RestoreV2:70*
use_locking(*
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/bn/beta*
validate_shape(*
_output_shapes
:@
Ò
save/Assign_71Assign!layerfilter1PL_edgefea_0/bn/gammasave/RestoreV2:71*
use_locking(*
T0*4
_class*
(&loc:@layerfilter1PL_edgefea_0/bn/gamma*
validate_shape(*
_output_shapes
:@
Ð
save/Assign_72Assign`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:72*
use_locking(*
T0*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
Ô
save/Assign_73Assignblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:73*
use_locking(*
T0*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
ì
save/Assign_74Assign layerfilter1PL_edgefea_0/weightssave/RestoreV2:74"/device:CPU:0*
use_locking(*
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
validate_shape(*'
_output_shapes
:@
ó
save/Assign_75Assign*layerfilter1PL_neib_att_conv_head_0/biasessave/RestoreV2:75"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter1PL_neib_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:
æ
save/Assign_76Assign+layerfilter1PL_neib_att_conv_head_0/bn/betasave/RestoreV2:76*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
è
save/Assign_77Assign,layerfilter1PL_neib_att_conv_head_0/bn/gammasave/RestoreV2:77*
use_locking(*
T0*?
_class5
31loc:@layerfilter1PL_neib_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
ý
save/Assign_78Assignvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:78*
use_locking(*
T0*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:

save/Assign_79Assignxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:79*
_output_shapes
:*
use_locking(*
T0*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(

save/Assign_80Assign+layerfilter1PL_neib_att_conv_head_0/weightssave/RestoreV2:80"/device:CPU:0*&
_output_shapes
:@*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
validate_shape(
â
save/Assign_81Assign)layerfilter1PL_newfea_conv_head_0/bn/betasave/RestoreV2:81*
use_locking(*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:@
ä
save/Assign_82Assign*layerfilter1PL_newfea_conv_head_0/bn/gammasave/RestoreV2:82*=
_class3
1/loc:@layerfilter1PL_newfea_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:@*
use_locking(*
T0
õ
save/Assign_83Assignrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:83*
use_locking(*
T0*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
ù
save/Assign_84Assigntlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:84*
use_locking(*
T0*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
þ
save/Assign_85Assign)layerfilter1PL_newfea_conv_head_0/weightssave/RestoreV2:85"/device:CPU:0*
use_locking(*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
validate_shape(*'
_output_shapes
:@
ó
save/Assign_86Assign*layerfilter1PL_self_att_conv_head_0/biasessave/RestoreV2:86"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter1PL_self_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:
æ
save/Assign_87Assign+layerfilter1PL_self_att_conv_head_0/bn/betasave/RestoreV2:87*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
è
save/Assign_88Assign,layerfilter1PL_self_att_conv_head_0/bn/gammasave/RestoreV2:88*
_output_shapes
:*
use_locking(*
T0*?
_class5
31loc:@layerfilter1PL_self_att_conv_head_0/bn/gamma*
validate_shape(
ý
save/Assign_89Assignvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:89*
use_locking(*
T0*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:

save/Assign_90Assignxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:90*
use_locking(*
T0*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:

save/Assign_91Assign+layerfilter1PL_self_att_conv_head_0/weightssave/RestoreV2:91"/device:CPU:0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
validate_shape(*&
_output_shapes
:@*
use_locking(*
T0
´
save/restore_all/NoOpNoOp^save/Assign_1^save/Assign_10^save/Assign_13^save/Assign_14^save/Assign_15^save/Assign_16^save/Assign_2^save/Assign_21^save/Assign_22^save/Assign_23^save/Assign_24^save/Assign_27^save/Assign_28^save/Assign_29^save/Assign_3^save/Assign_30^save/Assign_33^save/Assign_34^save/Assign_35^save/Assign_36^save/Assign_39^save/Assign_4^save/Assign_40^save/Assign_41^save/Assign_42^save/Assign_44^save/Assign_46^save/Assign_47^save/Assign_48^save/Assign_49^save/Assign_52^save/Assign_53^save/Assign_54^save/Assign_55^save/Assign_57^save/Assign_58^save/Assign_59^save/Assign_60^save/Assign_63^save/Assign_64^save/Assign_65^save/Assign_66^save/Assign_68^save/Assign_7^save/Assign_70^save/Assign_71^save/Assign_72^save/Assign_73^save/Assign_76^save/Assign_77^save/Assign_78^save/Assign_79^save/Assign_8^save/Assign_81^save/Assign_82^save/Assign_83^save/Assign_84^save/Assign_87^save/Assign_88^save/Assign_89^save/Assign_9^save/Assign_90
§
save/restore_all/NoOp_1NoOp^save/Assign^save/Assign_11^save/Assign_12^save/Assign_17^save/Assign_18^save/Assign_19^save/Assign_20^save/Assign_25^save/Assign_26^save/Assign_31^save/Assign_32^save/Assign_37^save/Assign_38^save/Assign_43^save/Assign_45^save/Assign_5^save/Assign_50^save/Assign_51^save/Assign_56^save/Assign_6^save/Assign_61^save/Assign_62^save/Assign_67^save/Assign_69^save/Assign_74^save/Assign_75^save/Assign_80^save/Assign_85^save/Assign_86^save/Assign_91"/device:CPU:0
J
save/restore_allNoOp^save/restore_all/NoOp^save/restore_all/NoOp_1
[
save_1/filename/inputConst*
_output_shapes
: *
valueB Bmodel*
dtype0
r
save_1/filenamePlaceholderWithDefaultsave_1/filename/input*
dtype0*
_output_shapes
: *
shape: 
i
save_1/ConstPlaceholderWithDefaultsave_1/filename*
_output_shapes
: *
shape: *
dtype0

save_1/StringJoin/inputs_1Const*<
value3B1 B+_temp_5eabc7e096284e9881e9e106b622bb3a/part*
dtype0*
_output_shapes
: 
{
save_1/StringJoin
StringJoinsave_1/Constsave_1/StringJoin/inputs_1*
	separator *
N*
_output_shapes
: 
S
save_1/num_shardsConst*
_output_shapes
: *
value	B :*
dtype0
m
save_1/ShardedFilename/shardConst"/device:CPU:0*
value	B : *
dtype0*
_output_shapes
: 

save_1/ShardedFilenameShardedFilenamesave_1/StringJoinsave_1/ShardedFilename/shardsave_1/num_shards"/device:CPU:0*
_output_shapes
: 

save_1/SaveV2/tensor_namesConst"/device:CPU:0*®
value¤B¡>B:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc1/bn/betaBPLfc1/bn/gammaB:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc2/bn/betaBPLfc2/bn/gammaB:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverageB<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageBaggPL/bn/betaBaggPL/bn/gammaBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01PL/bn/betaBgapnet01PL/bn/gammaBDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11PL/bn/betaBgapnet11PL/bn/gammaBDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter0PL/BiasAdd/biasesB layerfilter0PL_edgefea_0/bn/betaB!layerfilter0PL_edgefea_0/bn/gammaB`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_neib_att_conv_head_0/bn/betaB,layerfilter0PL_neib_att_conv_head_0/bn/gammaBvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0PL_newfea_conv_head_0/bn/betaB*layerfilter0PL_newfea_conv_head_0/bn/gammaBrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_self_att_conv_head_0/bn/betaB,layerfilter0PL_self_att_conv_head_0/bn/gammaBvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1PL/BiasAdd/biasesB layerfilter1PL_edgefea_0/bn/betaB!layerfilter1PL_edgefea_0/bn/gammaB`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_neib_att_conv_head_0/bn/betaB,layerfilter1PL_neib_att_conv_head_0/bn/gammaBvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1PL_newfea_conv_head_0/bn/betaB*layerfilter1PL_newfea_conv_head_0/bn/gammaBrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_self_att_conv_head_0/bn/betaB,layerfilter1PL_self_att_conv_head_0/bn/gammaBvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:>
ó
save_1/SaveV2/shape_and_slicesConst"/device:CPU:0*
valueB>B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:>
á
save_1/SaveV2SaveV2save_1/ShardedFilenamesave_1/SaveV2/tensor_namessave_1/SaveV2/shape_and_slices:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAveragePLfc1/bn/betaPLfc1/bn/gamma:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAveragePLfc2/bn/betaPLfc2/bn/gamma:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageaggPL/bn/betaaggPL/bn/gammagapnet00/bn/betagapnet00/bn/gamma@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet01PL/bn/betagapnet01PL/bn/gammaDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverageFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet10/bn/betagapnet10/bn/gamma@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet11PL/bn/betagapnet11PL/bn/gammaDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverageFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0PL/BiasAdd/biases layerfilter0PL_edgefea_0/bn/beta!layerfilter0PL_edgefea_0/bn/gamma`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter0PL_neib_att_conv_head_0/bn/beta,layerfilter0PL_neib_att_conv_head_0/bn/gammavlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter0PL_newfea_conv_head_0/bn/beta*layerfilter0PL_newfea_conv_head_0/bn/gammarlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter0PL_self_att_conv_head_0/bn/beta,layerfilter0PL_self_att_conv_head_0/bn/gammavlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1PL/BiasAdd/biases layerfilter1PL_edgefea_0/bn/beta!layerfilter1PL_edgefea_0/bn/gamma`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter1PL_neib_att_conv_head_0/bn/beta,layerfilter1PL_neib_att_conv_head_0/bn/gammavlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1PL_newfea_conv_head_0/bn/beta*layerfilter1PL_newfea_conv_head_0/bn/gammarlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage+layerfilter1PL_self_att_conv_head_0/bn/beta,layerfilter1PL_self_att_conv_head_0/bn/gammavlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragexlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage"/device:CPU:0*L
dtypesB
@2>
¨
save_1/control_dependencyIdentitysave_1/ShardedFilename^save_1/SaveV2"/device:CPU:0*)
_class
loc:@save_1/ShardedFilename*
_output_shapes
: *
T0
o
save_1/ShardedFilename_1/shardConst"/device:CPU:0*
value	B :*
dtype0*
_output_shapes
: 

save_1/ShardedFilename_1ShardedFilenamesave_1/StringJoinsave_1/ShardedFilename_1/shardsave_1/num_shards"/device:CPU:0*
_output_shapes
: 
½
save_1/SaveV2_1/tensor_namesConst"/device:CPU:0*Ý
valueÓBÐBPLfc1/biasesBPLfc1/weightsBPLfc2/biasesBPLfc2/weightsBaggPL/biasesBaggPL/weightsBfc3PL/biasesBfc3PL/weightsBgapnet00/biasesBgapnet00/weightsBgapnet01PL/biasesBgapnet01PL/weightsBgapnet10/biasesBgapnet10/weightsBgapnet11PL/biasesBgapnet11PL/weightsBlayerfilter0PL_edgefea_0/biasesB layerfilter0PL_edgefea_0/weightsB*layerfilter0PL_neib_att_conv_head_0/biasesB+layerfilter0PL_neib_att_conv_head_0/weightsB)layerfilter0PL_newfea_conv_head_0/weightsB*layerfilter0PL_self_att_conv_head_0/biasesB+layerfilter0PL_self_att_conv_head_0/weightsBlayerfilter1PL_edgefea_0/biasesB layerfilter1PL_edgefea_0/weightsB*layerfilter1PL_neib_att_conv_head_0/biasesB+layerfilter1PL_neib_att_conv_head_0/weightsB)layerfilter1PL_newfea_conv_head_0/weightsB*layerfilter1PL_self_att_conv_head_0/biasesB+layerfilter1PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
:
²
 save_1/SaveV2_1/shape_and_slicesConst"/device:CPU:0*O
valueFBDB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
ø
save_1/SaveV2_1SaveV2save_1/ShardedFilename_1save_1/SaveV2_1/tensor_names save_1/SaveV2_1/shape_and_slicesPLfc1/biasesPLfc1/weightsPLfc2/biasesPLfc2/weightsaggPL/biasesaggPL/weightsfc3PL/biasesfc3PL/weightsgapnet00/biasesgapnet00/weightsgapnet01PL/biasesgapnet01PL/weightsgapnet10/biasesgapnet10/weightsgapnet11PL/biasesgapnet11PL/weightslayerfilter0PL_edgefea_0/biases layerfilter0PL_edgefea_0/weights*layerfilter0PL_neib_att_conv_head_0/biases+layerfilter0PL_neib_att_conv_head_0/weights)layerfilter0PL_newfea_conv_head_0/weights*layerfilter0PL_self_att_conv_head_0/biases+layerfilter0PL_self_att_conv_head_0/weightslayerfilter1PL_edgefea_0/biases layerfilter1PL_edgefea_0/weights*layerfilter1PL_neib_att_conv_head_0/biases+layerfilter1PL_neib_att_conv_head_0/weights)layerfilter1PL_newfea_conv_head_0/weights*layerfilter1PL_self_att_conv_head_0/biases+layerfilter1PL_self_att_conv_head_0/weights"/device:CPU:0*,
dtypes"
 2
°
save_1/control_dependency_1Identitysave_1/ShardedFilename_1^save_1/SaveV2_1"/device:CPU:0*
T0*+
_class!
loc:@save_1/ShardedFilename_1*
_output_shapes
: 
ê
-save_1/MergeV2Checkpoints/checkpoint_prefixesPacksave_1/ShardedFilenamesave_1/ShardedFilename_1^save_1/control_dependency^save_1/control_dependency_1"/device:CPU:0*
T0*

axis *
N*
_output_shapes
:

save_1/MergeV2CheckpointsMergeV2Checkpoints-save_1/MergeV2Checkpoints/checkpoint_prefixessave_1/Const"/device:CPU:0*
delete_old_dirs(
¯
save_1/IdentityIdentitysave_1/Const^save_1/MergeV2Checkpoints^save_1/control_dependency^save_1/control_dependency_1"/device:CPU:0*
_output_shapes
: *
T0

save_1/RestoreV2/tensor_namesConst"/device:CPU:0*®
value¤B¡>B:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc1/bn/betaBPLfc1/bn/gammaB:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverageB<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverageBPLfc2/bn/betaBPLfc2/bn/gammaB:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverageB<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverageBaggPL/bn/betaBaggPL/bn/gammaBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01PL/bn/betaBgapnet01PL/bn/gammaBDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11PL/bn/betaBgapnet11PL/bn/gammaBDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverageBFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter0PL/BiasAdd/biasesB layerfilter0PL_edgefea_0/bn/betaB!layerfilter0PL_edgefea_0/bn/gammaB`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_neib_att_conv_head_0/bn/betaB,layerfilter0PL_neib_att_conv_head_0/bn/gammaBvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0PL_newfea_conv_head_0/bn/betaB*layerfilter0PL_newfea_conv_head_0/bn/gammaBrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter0PL_self_att_conv_head_0/bn/betaB,layerfilter0PL_self_att_conv_head_0/bn/gammaBvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1PL/BiasAdd/biasesB layerfilter1PL_edgefea_0/bn/betaB!layerfilter1PL_edgefea_0/bn/gammaB`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageBblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_neib_att_conv_head_0/bn/betaB,layerfilter1PL_neib_att_conv_head_0/bn/gammaBvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1PL_newfea_conv_head_0/bn/betaB*layerfilter1PL_newfea_conv_head_0/bn/gammaBrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB+layerfilter1PL_self_att_conv_head_0/bn/betaB,layerfilter1PL_self_att_conv_head_0/bn/gammaBvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:>
ö
!save_1/RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:>*
valueB>B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0
Û
save_1/RestoreV2	RestoreV2save_1/Constsave_1/RestoreV2/tensor_names!save_1/RestoreV2/shape_and_slices"/device:CPU:0*
_output_shapesû
ø::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*L
dtypesB
@2>

save_1/AssignAssign:PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2*
use_locking(*
T0*M
_classC
A?loc:@PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save_1/Assign_1Assign<PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:1*
use_locking(*
T0*O
_classE
CAloc:@PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
­
save_1/Assign_2AssignPLfc1/bn/betasave_1/RestoreV2:2*
_output_shapes	
:*
use_locking(*
T0* 
_class
loc:@PLfc1/bn/beta*
validate_shape(
¯
save_1/Assign_3AssignPLfc1/bn/gammasave_1/RestoreV2:3*!
_class
loc:@PLfc1/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

save_1/Assign_4Assign:PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:4*
use_locking(*
T0*M
_classC
A?loc:@PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save_1/Assign_5Assign<PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:5*
use_locking(*
T0*O
_classE
CAloc:@PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
­
save_1/Assign_6AssignPLfc2/bn/betasave_1/RestoreV2:6*
_output_shapes	
:*
use_locking(*
T0* 
_class
loc:@PLfc2/bn/beta*
validate_shape(
¯
save_1/Assign_7AssignPLfc2/bn/gammasave_1/RestoreV2:7*
use_locking(*
T0*!
_class
loc:@PLfc2/bn/gamma*
validate_shape(*
_output_shapes	
:

save_1/Assign_8Assign:aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:8*
use_locking(*
T0*M
_classC
A?loc:@aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save_1/Assign_9Assign<aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:9*O
_classE
CAloc:@aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
¯
save_1/Assign_10AssignaggPL/bn/betasave_1/RestoreV2:10*
use_locking(*
T0* 
_class
loc:@aggPL/bn/beta*
validate_shape(*
_output_shapes	
:
±
save_1/Assign_11AssignaggPL/bn/gammasave_1/RestoreV2:11*
use_locking(*
T0*!
_class
loc:@aggPL/bn/gamma*
validate_shape(*
_output_shapes	
:
´
save_1/Assign_12Assigngapnet00/bn/betasave_1/RestoreV2:12*
use_locking(*
T0*#
_class
loc:@gapnet00/bn/beta*
validate_shape(*
_output_shapes
:@
¶
save_1/Assign_13Assigngapnet00/bn/gammasave_1/RestoreV2:13*$
_class
loc:@gapnet00/bn/gamma*
validate_shape(*
_output_shapes
:@*
use_locking(*
T0

save_1/Assign_14Assign@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:14*
use_locking(*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@

save_1/Assign_15AssignBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:15*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@*
use_locking(*
T0
¹
save_1/Assign_16Assigngapnet01PL/bn/betasave_1/RestoreV2:16*
_output_shapes	
:*
use_locking(*
T0*%
_class
loc:@gapnet01PL/bn/beta*
validate_shape(
»
save_1/Assign_17Assigngapnet01PL/bn/gammasave_1/RestoreV2:17*&
_class
loc:@gapnet01PL/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

save_1/Assign_18AssignDgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:18*
use_locking(*
T0*W
_classM
KIloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
¡
save_1/Assign_19AssignFgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:19*Y
_classO
MKloc:@gapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
µ
save_1/Assign_20Assigngapnet10/bn/betasave_1/RestoreV2:20*
use_locking(*
T0*#
_class
loc:@gapnet10/bn/beta*
validate_shape(*
_output_shapes	
:
·
save_1/Assign_21Assigngapnet10/bn/gammasave_1/RestoreV2:21*
use_locking(*
T0*$
_class
loc:@gapnet10/bn/gamma*
validate_shape(*
_output_shapes	
:

save_1/Assign_22Assign@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:22*
use_locking(*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:

save_1/Assign_23AssignBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:23*
use_locking(*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
¹
save_1/Assign_24Assigngapnet11PL/bn/betasave_1/RestoreV2:24*%
_class
loc:@gapnet11PL/bn/beta*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
»
save_1/Assign_25Assigngapnet11PL/bn/gammasave_1/RestoreV2:25*&
_class
loc:@gapnet11PL/bn/gamma*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0

save_1/Assign_26AssignDgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:26*
use_locking(*
T0*W
_classM
KIloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
¡
save_1/Assign_27AssignFgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:27*
use_locking(*
T0*Y
_classO
MKloc:@gapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes	
:
Î
save_1/Assign_28Assignlayerfilter0PL/BiasAdd/biasessave_1/RestoreV2:28*0
_class&
$"loc:@layerfilter0PL/BiasAdd/biases*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
Ô
save_1/Assign_29Assign layerfilter0PL_edgefea_0/bn/betasave_1/RestoreV2:29*3
_class)
'%loc:@layerfilter0PL_edgefea_0/bn/beta*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
Ö
save_1/Assign_30Assign!layerfilter0PL_edgefea_0/bn/gammasave_1/RestoreV2:30*
use_locking(*
T0*4
_class*
(&loc:@layerfilter0PL_edgefea_0/bn/gamma*
validate_shape(*
_output_shapes
:
Ô
save_1/Assign_31Assign`layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:31*s
_classi
geloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
Ø
save_1/Assign_32Assignblayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:32*
use_locking(*
T0*u
_classk
igloc:@layerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ê
save_1/Assign_33Assign+layerfilter0PL_neib_att_conv_head_0/bn/betasave_1/RestoreV2:33*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
ì
save_1/Assign_34Assign,layerfilter0PL_neib_att_conv_head_0/bn/gammasave_1/RestoreV2:34*?
_class5
31loc:@layerfilter0PL_neib_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:*
use_locking(*
T0

save_1/Assign_35Assignvlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:35*
_output_shapes
:*
use_locking(*
T0*
_class
}{loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(

save_1/Assign_36Assignxlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:36*
_output_shapes
:*
use_locking(*
T0*
_class
}loc:@layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
æ
save_1/Assign_37Assign)layerfilter0PL_newfea_conv_head_0/bn/betasave_1/RestoreV2:37*
use_locking(*
T0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
è
save_1/Assign_38Assign*layerfilter0PL_newfea_conv_head_0/bn/gammasave_1/RestoreV2:38*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_newfea_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:
ù
save_1/Assign_39Assignrlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:39*
use_locking(*
T0*
_class{
ywloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ý
save_1/Assign_40Assigntlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:40*
use_locking(*
T0*
_class}
{yloc:@layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ê
save_1/Assign_41Assign+layerfilter0PL_self_att_conv_head_0/bn/betasave_1/RestoreV2:41*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:
ì
save_1/Assign_42Assign,layerfilter0PL_self_att_conv_head_0/bn/gammasave_1/RestoreV2:42*?
_class5
31loc:@layerfilter0PL_self_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:*
use_locking(*
T0

save_1/Assign_43Assignvlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:43*
use_locking(*
T0*
_class
}{loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:

save_1/Assign_44Assignxlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:44*
use_locking(*
T0*
_class
}loc:@layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
Î
save_1/Assign_45Assignlayerfilter1PL/BiasAdd/biasessave_1/RestoreV2:45*
use_locking(*
T0*0
_class&
$"loc:@layerfilter1PL/BiasAdd/biases*
validate_shape(*
_output_shapes
:@
Ô
save_1/Assign_46Assign layerfilter1PL_edgefea_0/bn/betasave_1/RestoreV2:46*
use_locking(*
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/bn/beta*
validate_shape(*
_output_shapes
:@
Ö
save_1/Assign_47Assign!layerfilter1PL_edgefea_0/bn/gammasave_1/RestoreV2:47*
use_locking(*
T0*4
_class*
(&loc:@layerfilter1PL_edgefea_0/bn/gamma*
validate_shape(*
_output_shapes
:@
Ô
save_1/Assign_48Assign`layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:48*
_output_shapes
:@*
use_locking(*
T0*s
_classi
geloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(
Ø
save_1/Assign_49Assignblayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:49*
use_locking(*
T0*u
_classk
igloc:@layerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
ê
save_1/Assign_50Assign+layerfilter1PL_neib_att_conv_head_0/bn/betasave_1/RestoreV2:50*
_output_shapes
:*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/bn/beta*
validate_shape(
ì
save_1/Assign_51Assign,layerfilter1PL_neib_att_conv_head_0/bn/gammasave_1/RestoreV2:51*?
_class5
31loc:@layerfilter1PL_neib_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:*
use_locking(*
T0

save_1/Assign_52Assignvlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:52*
_output_shapes
:*
use_locking(*
T0*
_class
}{loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(

save_1/Assign_53Assignxlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:53*
_output_shapes
:*
use_locking(*
T0*
_class
}loc:@layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
æ
save_1/Assign_54Assign)layerfilter1PL_newfea_conv_head_0/bn/betasave_1/RestoreV2:54*
use_locking(*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/bn/beta*
validate_shape(*
_output_shapes
:@
è
save_1/Assign_55Assign*layerfilter1PL_newfea_conv_head_0/bn/gammasave_1/RestoreV2:55*
_output_shapes
:@*
use_locking(*
T0*=
_class3
1/loc:@layerfilter1PL_newfea_conv_head_0/bn/gamma*
validate_shape(
ù
save_1/Assign_56Assignrlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:56*
use_locking(*
T0*
_class{
ywloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:@
ý
save_1/Assign_57Assigntlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:57*
_output_shapes
:@*
use_locking(*
T0*
_class}
{yloc:@layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(
ê
save_1/Assign_58Assign+layerfilter1PL_self_att_conv_head_0/bn/betasave_1/RestoreV2:58*
_output_shapes
:*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/bn/beta*
validate_shape(
ì
save_1/Assign_59Assign,layerfilter1PL_self_att_conv_head_0/bn/gammasave_1/RestoreV2:59*
use_locking(*
T0*?
_class5
31loc:@layerfilter1PL_self_att_conv_head_0/bn/gamma*
validate_shape(*
_output_shapes
:

save_1/Assign_60Assignvlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:60*
_class
}{loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:*
use_locking(*
T0

save_1/Assign_61Assignxlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:61*
use_locking(*
T0*
_class
}loc:@layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
validate_shape(*
_output_shapes
:
ª	
save_1/restore_shardNoOp^save_1/Assign^save_1/Assign_1^save_1/Assign_10^save_1/Assign_11^save_1/Assign_12^save_1/Assign_13^save_1/Assign_14^save_1/Assign_15^save_1/Assign_16^save_1/Assign_17^save_1/Assign_18^save_1/Assign_19^save_1/Assign_2^save_1/Assign_20^save_1/Assign_21^save_1/Assign_22^save_1/Assign_23^save_1/Assign_24^save_1/Assign_25^save_1/Assign_26^save_1/Assign_27^save_1/Assign_28^save_1/Assign_29^save_1/Assign_3^save_1/Assign_30^save_1/Assign_31^save_1/Assign_32^save_1/Assign_33^save_1/Assign_34^save_1/Assign_35^save_1/Assign_36^save_1/Assign_37^save_1/Assign_38^save_1/Assign_39^save_1/Assign_4^save_1/Assign_40^save_1/Assign_41^save_1/Assign_42^save_1/Assign_43^save_1/Assign_44^save_1/Assign_45^save_1/Assign_46^save_1/Assign_47^save_1/Assign_48^save_1/Assign_49^save_1/Assign_5^save_1/Assign_50^save_1/Assign_51^save_1/Assign_52^save_1/Assign_53^save_1/Assign_54^save_1/Assign_55^save_1/Assign_56^save_1/Assign_57^save_1/Assign_58^save_1/Assign_59^save_1/Assign_6^save_1/Assign_60^save_1/Assign_61^save_1/Assign_7^save_1/Assign_8^save_1/Assign_9
À
save_1/RestoreV2_1/tensor_namesConst"/device:CPU:0*Ý
valueÓBÐBPLfc1/biasesBPLfc1/weightsBPLfc2/biasesBPLfc2/weightsBaggPL/biasesBaggPL/weightsBfc3PL/biasesBfc3PL/weightsBgapnet00/biasesBgapnet00/weightsBgapnet01PL/biasesBgapnet01PL/weightsBgapnet10/biasesBgapnet10/weightsBgapnet11PL/biasesBgapnet11PL/weightsBlayerfilter0PL_edgefea_0/biasesB layerfilter0PL_edgefea_0/weightsB*layerfilter0PL_neib_att_conv_head_0/biasesB+layerfilter0PL_neib_att_conv_head_0/weightsB)layerfilter0PL_newfea_conv_head_0/weightsB*layerfilter0PL_self_att_conv_head_0/biasesB+layerfilter0PL_self_att_conv_head_0/weightsBlayerfilter1PL_edgefea_0/biasesB layerfilter1PL_edgefea_0/weightsB*layerfilter1PL_neib_att_conv_head_0/biasesB+layerfilter1PL_neib_att_conv_head_0/weightsB)layerfilter1PL_newfea_conv_head_0/weightsB*layerfilter1PL_self_att_conv_head_0/biasesB+layerfilter1PL_self_att_conv_head_0/weights*
dtype0*
_output_shapes
:
µ
#save_1/RestoreV2_1/shape_and_slicesConst"/device:CPU:0*O
valueFBDB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:
¿
save_1/RestoreV2_1	RestoreV2save_1/Constsave_1/RestoreV2_1/tensor_names#save_1/RestoreV2_1/shape_and_slices"/device:CPU:0*
_output_shapesz
x::::::::::::::::::::::::::::::*,
dtypes"
 2
»
save_1/Assign_62AssignPLfc1/biasessave_1/RestoreV2_1"/device:CPU:0*
use_locking(*
T0*
_class
loc:@PLfc1/biases*
validate_shape(*
_output_shapes	
:
Ä
save_1/Assign_63AssignPLfc1/weightssave_1/RestoreV2_1:1"/device:CPU:0* 
_output_shapes
:
*
use_locking(*
T0* 
_class
loc:@PLfc1/weights*
validate_shape(
½
save_1/Assign_64AssignPLfc2/biasessave_1/RestoreV2_1:2"/device:CPU:0*
use_locking(*
T0*
_class
loc:@PLfc2/biases*
validate_shape(*
_output_shapes	
:
Ä
save_1/Assign_65AssignPLfc2/weightssave_1/RestoreV2_1:3"/device:CPU:0* 
_class
loc:@PLfc2/weights*
validate_shape(* 
_output_shapes
:
*
use_locking(*
T0
½
save_1/Assign_66AssignaggPL/biasessave_1/RestoreV2_1:4"/device:CPU:0*
_class
loc:@aggPL/biases*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
Ì
save_1/Assign_67AssignaggPL/weightssave_1/RestoreV2_1:5"/device:CPU:0*
use_locking(*
T0* 
_class
loc:@aggPL/weights*
validate_shape(*(
_output_shapes
:
¼
save_1/Assign_68Assignfc3PL/biasessave_1/RestoreV2_1:6"/device:CPU:0*
_class
loc:@fc3PL/biases*
validate_shape(*
_output_shapes
:*
use_locking(*
T0
Ã
save_1/Assign_69Assignfc3PL/weightssave_1/RestoreV2_1:7"/device:CPU:0*
use_locking(*
T0* 
_class
loc:@fc3PL/weights*
validate_shape(*
_output_shapes
:	
Â
save_1/Assign_70Assigngapnet00/biasessave_1/RestoreV2_1:8"/device:CPU:0*
_output_shapes
:@*
use_locking(*
T0*"
_class
loc:@gapnet00/biases*
validate_shape(
Ð
save_1/Assign_71Assigngapnet00/weightssave_1/RestoreV2_1:9"/device:CPU:0*&
_output_shapes
:@*
use_locking(*
T0*#
_class
loc:@gapnet00/weights*
validate_shape(
È
save_1/Assign_72Assigngapnet01PL/biasessave_1/RestoreV2_1:10"/device:CPU:0*$
_class
loc:@gapnet01PL/biases*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
Ö
save_1/Assign_73Assigngapnet01PL/weightssave_1/RestoreV2_1:11"/device:CPU:0*%
_class
loc:@gapnet01PL/weights*
validate_shape(*'
_output_shapes
:@*
use_locking(*
T0
Ä
save_1/Assign_74Assigngapnet10/biasessave_1/RestoreV2_1:12"/device:CPU:0*
use_locking(*
T0*"
_class
loc:@gapnet10/biases*
validate_shape(*
_output_shapes	
:
Ò
save_1/Assign_75Assigngapnet10/weightssave_1/RestoreV2_1:13"/device:CPU:0*
use_locking(*
T0*#
_class
loc:@gapnet10/weights*
validate_shape(*'
_output_shapes
:M
È
save_1/Assign_76Assigngapnet11PL/biasessave_1/RestoreV2_1:14"/device:CPU:0*$
_class
loc:@gapnet11PL/biases*
validate_shape(*
_output_shapes	
:*
use_locking(*
T0
×
save_1/Assign_77Assigngapnet11PL/weightssave_1/RestoreV2_1:15"/device:CPU:0*(
_output_shapes
:*
use_locking(*
T0*%
_class
loc:@gapnet11PL/weights*
validate_shape(
ã
save_1/Assign_78Assignlayerfilter0PL_edgefea_0/biasessave_1/RestoreV2_1:16"/device:CPU:0*
_output_shapes
:*
use_locking(*
T0*2
_class(
&$loc:@layerfilter0PL_edgefea_0/biases*
validate_shape(
ñ
save_1/Assign_79Assign layerfilter0PL_edgefea_0/weightssave_1/RestoreV2_1:17"/device:CPU:0*3
_class)
'%loc:@layerfilter0PL_edgefea_0/weights*
validate_shape(*&
_output_shapes
:*
use_locking(*
T0
ù
save_1/Assign_80Assign*layerfilter0PL_neib_att_conv_head_0/biasessave_1/RestoreV2_1:18"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_neib_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:

save_1/Assign_81Assign+layerfilter0PL_neib_att_conv_head_0/weightssave_1/RestoreV2_1:19"/device:CPU:0*&
_output_shapes
:*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_neib_att_conv_head_0/weights*
validate_shape(

save_1/Assign_82Assign)layerfilter0PL_newfea_conv_head_0/weightssave_1/RestoreV2_1:20"/device:CPU:0*<
_class2
0.loc:@layerfilter0PL_newfea_conv_head_0/weights*
validate_shape(*&
_output_shapes
:*
use_locking(*
T0
ù
save_1/Assign_83Assign*layerfilter0PL_self_att_conv_head_0/biasessave_1/RestoreV2_1:21"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter0PL_self_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:

save_1/Assign_84Assign+layerfilter0PL_self_att_conv_head_0/weightssave_1/RestoreV2_1:22"/device:CPU:0*
use_locking(*
T0*>
_class4
20loc:@layerfilter0PL_self_att_conv_head_0/weights*
validate_shape(*&
_output_shapes
:
ã
save_1/Assign_85Assignlayerfilter1PL_edgefea_0/biasessave_1/RestoreV2_1:23"/device:CPU:0*2
_class(
&$loc:@layerfilter1PL_edgefea_0/biases*
validate_shape(*
_output_shapes
:@*
use_locking(*
T0
ò
save_1/Assign_86Assign layerfilter1PL_edgefea_0/weightssave_1/RestoreV2_1:24"/device:CPU:0*
use_locking(*
T0*3
_class)
'%loc:@layerfilter1PL_edgefea_0/weights*
validate_shape(*'
_output_shapes
:@
ù
save_1/Assign_87Assign*layerfilter1PL_neib_att_conv_head_0/biasessave_1/RestoreV2_1:25"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter1PL_neib_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:

save_1/Assign_88Assign+layerfilter1PL_neib_att_conv_head_0/weightssave_1/RestoreV2_1:26"/device:CPU:0*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_neib_att_conv_head_0/weights*
validate_shape(*&
_output_shapes
:@

save_1/Assign_89Assign)layerfilter1PL_newfea_conv_head_0/weightssave_1/RestoreV2_1:27"/device:CPU:0*
use_locking(*
T0*<
_class2
0.loc:@layerfilter1PL_newfea_conv_head_0/weights*
validate_shape(*'
_output_shapes
:@
ù
save_1/Assign_90Assign*layerfilter1PL_self_att_conv_head_0/biasessave_1/RestoreV2_1:28"/device:CPU:0*
use_locking(*
T0*=
_class3
1/loc:@layerfilter1PL_self_att_conv_head_0/biases*
validate_shape(*
_output_shapes
:

save_1/Assign_91Assign+layerfilter1PL_self_att_conv_head_0/weightssave_1/RestoreV2_1:29"/device:CPU:0*&
_output_shapes
:@*
use_locking(*
T0*>
_class4
20loc:@layerfilter1PL_self_att_conv_head_0/weights*
validate_shape(
ç
save_1/restore_shard_1NoOp^save_1/Assign_62^save_1/Assign_63^save_1/Assign_64^save_1/Assign_65^save_1/Assign_66^save_1/Assign_67^save_1/Assign_68^save_1/Assign_69^save_1/Assign_70^save_1/Assign_71^save_1/Assign_72^save_1/Assign_73^save_1/Assign_74^save_1/Assign_75^save_1/Assign_76^save_1/Assign_77^save_1/Assign_78^save_1/Assign_79^save_1/Assign_80^save_1/Assign_81^save_1/Assign_82^save_1/Assign_83^save_1/Assign_84^save_1/Assign_85^save_1/Assign_86^save_1/Assign_87^save_1/Assign_88^save_1/Assign_89^save_1/Assign_90^save_1/Assign_91"/device:CPU:0
6
save_1/restore_all/NoOpNoOp^save_1/restore_shard
I
save_1/restore_all/NoOp_1NoOp^save_1/restore_shard_1"/device:CPU:0
P
save_1/restore_allNoOp^save_1/restore_all/NoOp^save_1/restore_all/NoOp_1"&B
save_1/Const:0save_1/Identity:0save_1/restore_all (5 @F8"Ä
losses¹
¶
/layerfilter0PL_newfea_conv_head_0/weight_loss:0
&layerfilter0PL_edgefea_0/weight_loss:0
1layerfilter0PL_self_att_conv_head_0/weight_loss:0
1layerfilter0PL_neib_att_conv_head_0/weight_loss:0
gapnet00/weight_loss:0
gapnet01PL/weight_loss:0
/layerfilter1PL_newfea_conv_head_0/weight_loss:0
&layerfilter1PL_edgefea_0/weight_loss:0
1layerfilter1PL_self_att_conv_head_0/weight_loss:0
1layerfilter1PL_neib_att_conv_head_0/weight_loss:0
gapnet10/weight_loss:0
gapnet11PL/weight_loss:0
aggPL/weight_loss:0
PLfc1/weight_loss:0
PLfc2/weight_loss:0
fc3PL/weight_loss:0"á
model_variablesÍÊ
¢
layerfilter0PL/BiasAdd/biases:0$layerfilter0PL/BiasAdd/biases/Assign$layerfilter0PL/BiasAdd/biases/read:021layerfilter0PL/BiasAdd/biases/Initializer/zeros:08
¢
layerfilter1PL/BiasAdd/biases:0$layerfilter1PL/BiasAdd/biases/Assign$layerfilter1PL/BiasAdd/biases/read:021layerfilter1PL/BiasAdd/biases/Initializer/zeros:08"J
trainable_variablesJýI
Û
+layerfilter0PL_newfea_conv_head_0/weights:00layerfilter0PL_newfea_conv_head_0/weights/Assign0layerfilter0PL_newfea_conv_head_0/weights/read:02Flayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform:08
Á
+layerfilter0PL_newfea_conv_head_0/bn/beta:00layerfilter0PL_newfea_conv_head_0/bn/beta/Assign0layerfilter0PL_newfea_conv_head_0/bn/beta/read:02,layerfilter0PL_newfea_conv_head_0/bn/Const:08
Æ
,layerfilter0PL_newfea_conv_head_0/bn/gamma:01layerfilter0PL_newfea_conv_head_0/bn/gamma/Assign1layerfilter0PL_newfea_conv_head_0/bn/gamma/read:02.layerfilter0PL_newfea_conv_head_0/bn/Const_1:08
·
"layerfilter0PL_edgefea_0/weights:0'layerfilter0PL_edgefea_0/weights/Assign'layerfilter0PL_edgefea_0/weights/read:02=layerfilter0PL_edgefea_0/weights/Initializer/random_uniform:08
ª
!layerfilter0PL_edgefea_0/biases:0&layerfilter0PL_edgefea_0/biases/Assign&layerfilter0PL_edgefea_0/biases/read:023layerfilter0PL_edgefea_0/biases/Initializer/Const:08

"layerfilter0PL_edgefea_0/bn/beta:0'layerfilter0PL_edgefea_0/bn/beta/Assign'layerfilter0PL_edgefea_0/bn/beta/read:02#layerfilter0PL_edgefea_0/bn/Const:08
¢
#layerfilter0PL_edgefea_0/bn/gamma:0(layerfilter0PL_edgefea_0/bn/gamma/Assign(layerfilter0PL_edgefea_0/bn/gamma/read:02%layerfilter0PL_edgefea_0/bn/Const_1:08
ã
-layerfilter0PL_self_att_conv_head_0/weights:02layerfilter0PL_self_att_conv_head_0/weights/Assign2layerfilter0PL_self_att_conv_head_0/weights/read:02Hlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter0PL_self_att_conv_head_0/biases:01layerfilter0PL_self_att_conv_head_0/biases/Assign1layerfilter0PL_self_att_conv_head_0/biases/read:02>layerfilter0PL_self_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter0PL_self_att_conv_head_0/bn/beta:02layerfilter0PL_self_att_conv_head_0/bn/beta/Assign2layerfilter0PL_self_att_conv_head_0/bn/beta/read:02.layerfilter0PL_self_att_conv_head_0/bn/Const:08
Î
.layerfilter0PL_self_att_conv_head_0/bn/gamma:03layerfilter0PL_self_att_conv_head_0/bn/gamma/Assign3layerfilter0PL_self_att_conv_head_0/bn/gamma/read:020layerfilter0PL_self_att_conv_head_0/bn/Const_1:08
ã
-layerfilter0PL_neib_att_conv_head_0/weights:02layerfilter0PL_neib_att_conv_head_0/weights/Assign2layerfilter0PL_neib_att_conv_head_0/weights/read:02Hlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter0PL_neib_att_conv_head_0/biases:01layerfilter0PL_neib_att_conv_head_0/biases/Assign1layerfilter0PL_neib_att_conv_head_0/biases/read:02>layerfilter0PL_neib_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter0PL_neib_att_conv_head_0/bn/beta:02layerfilter0PL_neib_att_conv_head_0/bn/beta/Assign2layerfilter0PL_neib_att_conv_head_0/bn/beta/read:02.layerfilter0PL_neib_att_conv_head_0/bn/Const:08
Î
.layerfilter0PL_neib_att_conv_head_0/bn/gamma:03layerfilter0PL_neib_att_conv_head_0/bn/gamma/Assign3layerfilter0PL_neib_att_conv_head_0/bn/gamma/read:020layerfilter0PL_neib_att_conv_head_0/bn/Const_1:08
¢
layerfilter0PL/BiasAdd/biases:0$layerfilter0PL/BiasAdd/biases/Assign$layerfilter0PL/BiasAdd/biases/read:021layerfilter0PL/BiasAdd/biases/Initializer/zeros:08
w
gapnet00/weights:0gapnet00/weights/Assigngapnet00/weights/read:02-gapnet00/weights/Initializer/random_uniform:08
j
gapnet00/biases:0gapnet00/biases/Assigngapnet00/biases/read:02#gapnet00/biases/Initializer/Const:08
]
gapnet00/bn/beta:0gapnet00/bn/beta/Assigngapnet00/bn/beta/read:02gapnet00/bn/Const:08
b
gapnet00/bn/gamma:0gapnet00/bn/gamma/Assigngapnet00/bn/gamma/read:02gapnet00/bn/Const_1:08

gapnet01PL/weights:0gapnet01PL/weights/Assigngapnet01PL/weights/read:02/gapnet01PL/weights/Initializer/random_uniform:08
r
gapnet01PL/biases:0gapnet01PL/biases/Assigngapnet01PL/biases/read:02%gapnet01PL/biases/Initializer/Const:08
e
gapnet01PL/bn/beta:0gapnet01PL/bn/beta/Assigngapnet01PL/bn/beta/read:02gapnet01PL/bn/Const:08
j
gapnet01PL/bn/gamma:0gapnet01PL/bn/gamma/Assigngapnet01PL/bn/gamma/read:02gapnet01PL/bn/Const_1:08
Û
+layerfilter1PL_newfea_conv_head_0/weights:00layerfilter1PL_newfea_conv_head_0/weights/Assign0layerfilter1PL_newfea_conv_head_0/weights/read:02Flayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform:08
Á
+layerfilter1PL_newfea_conv_head_0/bn/beta:00layerfilter1PL_newfea_conv_head_0/bn/beta/Assign0layerfilter1PL_newfea_conv_head_0/bn/beta/read:02,layerfilter1PL_newfea_conv_head_0/bn/Const:08
Æ
,layerfilter1PL_newfea_conv_head_0/bn/gamma:01layerfilter1PL_newfea_conv_head_0/bn/gamma/Assign1layerfilter1PL_newfea_conv_head_0/bn/gamma/read:02.layerfilter1PL_newfea_conv_head_0/bn/Const_1:08
·
"layerfilter1PL_edgefea_0/weights:0'layerfilter1PL_edgefea_0/weights/Assign'layerfilter1PL_edgefea_0/weights/read:02=layerfilter1PL_edgefea_0/weights/Initializer/random_uniform:08
ª
!layerfilter1PL_edgefea_0/biases:0&layerfilter1PL_edgefea_0/biases/Assign&layerfilter1PL_edgefea_0/biases/read:023layerfilter1PL_edgefea_0/biases/Initializer/Const:08

"layerfilter1PL_edgefea_0/bn/beta:0'layerfilter1PL_edgefea_0/bn/beta/Assign'layerfilter1PL_edgefea_0/bn/beta/read:02#layerfilter1PL_edgefea_0/bn/Const:08
¢
#layerfilter1PL_edgefea_0/bn/gamma:0(layerfilter1PL_edgefea_0/bn/gamma/Assign(layerfilter1PL_edgefea_0/bn/gamma/read:02%layerfilter1PL_edgefea_0/bn/Const_1:08
ã
-layerfilter1PL_self_att_conv_head_0/weights:02layerfilter1PL_self_att_conv_head_0/weights/Assign2layerfilter1PL_self_att_conv_head_0/weights/read:02Hlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter1PL_self_att_conv_head_0/biases:01layerfilter1PL_self_att_conv_head_0/biases/Assign1layerfilter1PL_self_att_conv_head_0/biases/read:02>layerfilter1PL_self_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter1PL_self_att_conv_head_0/bn/beta:02layerfilter1PL_self_att_conv_head_0/bn/beta/Assign2layerfilter1PL_self_att_conv_head_0/bn/beta/read:02.layerfilter1PL_self_att_conv_head_0/bn/Const:08
Î
.layerfilter1PL_self_att_conv_head_0/bn/gamma:03layerfilter1PL_self_att_conv_head_0/bn/gamma/Assign3layerfilter1PL_self_att_conv_head_0/bn/gamma/read:020layerfilter1PL_self_att_conv_head_0/bn/Const_1:08
ã
-layerfilter1PL_neib_att_conv_head_0/weights:02layerfilter1PL_neib_att_conv_head_0/weights/Assign2layerfilter1PL_neib_att_conv_head_0/weights/read:02Hlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter1PL_neib_att_conv_head_0/biases:01layerfilter1PL_neib_att_conv_head_0/biases/Assign1layerfilter1PL_neib_att_conv_head_0/biases/read:02>layerfilter1PL_neib_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter1PL_neib_att_conv_head_0/bn/beta:02layerfilter1PL_neib_att_conv_head_0/bn/beta/Assign2layerfilter1PL_neib_att_conv_head_0/bn/beta/read:02.layerfilter1PL_neib_att_conv_head_0/bn/Const:08
Î
.layerfilter1PL_neib_att_conv_head_0/bn/gamma:03layerfilter1PL_neib_att_conv_head_0/bn/gamma/Assign3layerfilter1PL_neib_att_conv_head_0/bn/gamma/read:020layerfilter1PL_neib_att_conv_head_0/bn/Const_1:08
¢
layerfilter1PL/BiasAdd/biases:0$layerfilter1PL/BiasAdd/biases/Assign$layerfilter1PL/BiasAdd/biases/read:021layerfilter1PL/BiasAdd/biases/Initializer/zeros:08
w
gapnet10/weights:0gapnet10/weights/Assigngapnet10/weights/read:02-gapnet10/weights/Initializer/random_uniform:08
j
gapnet10/biases:0gapnet10/biases/Assigngapnet10/biases/read:02#gapnet10/biases/Initializer/Const:08
]
gapnet10/bn/beta:0gapnet10/bn/beta/Assigngapnet10/bn/beta/read:02gapnet10/bn/Const:08
b
gapnet10/bn/gamma:0gapnet10/bn/gamma/Assigngapnet10/bn/gamma/read:02gapnet10/bn/Const_1:08

gapnet11PL/weights:0gapnet11PL/weights/Assigngapnet11PL/weights/read:02/gapnet11PL/weights/Initializer/random_uniform:08
r
gapnet11PL/biases:0gapnet11PL/biases/Assigngapnet11PL/biases/read:02%gapnet11PL/biases/Initializer/Const:08
e
gapnet11PL/bn/beta:0gapnet11PL/bn/beta/Assigngapnet11PL/bn/beta/read:02gapnet11PL/bn/Const:08
j
gapnet11PL/bn/gamma:0gapnet11PL/bn/gamma/Assigngapnet11PL/bn/gamma/read:02gapnet11PL/bn/Const_1:08
k
aggPL/weights:0aggPL/weights/AssignaggPL/weights/read:02*aggPL/weights/Initializer/random_uniform:08
^
aggPL/biases:0aggPL/biases/AssignaggPL/biases/read:02 aggPL/biases/Initializer/Const:08
Q
aggPL/bn/beta:0aggPL/bn/beta/AssignaggPL/bn/beta/read:02aggPL/bn/Const:08
V
aggPL/bn/gamma:0aggPL/bn/gamma/AssignaggPL/bn/gamma/read:02aggPL/bn/Const_1:08
k
PLfc1/weights:0PLfc1/weights/AssignPLfc1/weights/read:02*PLfc1/weights/Initializer/random_uniform:08
^
PLfc1/biases:0PLfc1/biases/AssignPLfc1/biases/read:02 PLfc1/biases/Initializer/Const:08
Q
PLfc1/bn/beta:0PLfc1/bn/beta/AssignPLfc1/bn/beta/read:02PLfc1/bn/Const:08
V
PLfc1/bn/gamma:0PLfc1/bn/gamma/AssignPLfc1/bn/gamma/read:02PLfc1/bn/Const_1:08
k
PLfc2/weights:0PLfc2/weights/AssignPLfc2/weights/read:02*PLfc2/weights/Initializer/random_uniform:08
^
PLfc2/biases:0PLfc2/biases/AssignPLfc2/biases/read:02 PLfc2/biases/Initializer/Const:08
Q
PLfc2/bn/beta:0PLfc2/bn/beta/AssignPLfc2/bn/beta/read:02PLfc2/bn/Const:08
V
PLfc2/bn/gamma:0PLfc2/bn/gamma/AssignPLfc2/bn/gamma/read:02PLfc2/bn/Const_1:08
k
fc3PL/weights:0fc3PL/weights/Assignfc3PL/weights/read:02*fc3PL/weights/Initializer/random_uniform:08
^
fc3PL/biases:0fc3PL/biases/Assignfc3PL/biases/read:02 fc3PL/biases/Initializer/Const:08"·ü
cond_context¥ü¡ü
ø
3layerfilter0PL_newfea_conv_head_0/bn/cond/cond_text3layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:04layerfilter0PL_newfea_conv_head_0/bn/cond/switch_t:0 *Ó
[layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter0PL_newfea_conv_head_0/bn/cond/control_dependency:0
3layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:0
4layerfilter0PL_newfea_conv_head_0/bn/cond/switch_t:0
ylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze:0
8layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1:0
6layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze:0clayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1j
3layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:03layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:0¡
8layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1:0elayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1×
vlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1â
{layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Ó
tlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Þ
ylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter0PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
¾
5layerfilter0PL_newfea_conv_head_0/bn/cond/cond_text_13layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:04layerfilter0PL_newfea_conv_head_0/bn/cond/switch_f:0*
@layerfilter0PL_newfea_conv_head_0/bn/cond/control_dependency_1:0
3layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:0
4layerfilter0PL_newfea_conv_head_0/bn/cond/switch_f:0j
3layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:03layerfilter0PL_newfea_conv_head_0/bn/cond/pred_id:0
ä
5layerfilter0PL_newfea_conv_head_0/bn/cond_1/cond_text5layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:06layerfilter0PL_newfea_conv_head_0/bn/cond_1/switch_t:0 *¹
=layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity/Switch:1
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity:0
?layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1
8layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity_1:0
5layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:0
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/switch_t:0
6layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze:0
8layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1:0w
6layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze:0=layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity/Switch:1{
8layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1:0?layerfilter0PL_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1n
5layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:05layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:0
Ð	
7layerfilter0PL_newfea_conv_head_0/bn/cond_1/cond_text_15layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:06layerfilter0PL_newfea_conv_head_0/bn/cond_1/switch_f:0*¥
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_1:0
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_1:1
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_2:0
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_2:1
5layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:0
6layerfilter0PL_newfea_conv_head_0/bn/cond_1/switch_f:0
ylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_1:0n
5layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:05layerfilter0PL_newfea_conv_head_0/bn/cond_1/pred_id:0µ
{layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter0PL_newfea_conv_head_0/bn/cond_1/Switch_2:0
­
*layerfilter0PL_edgefea_0/bn/cond/cond_text*layerfilter0PL_edgefea_0/bn/cond/pred_id:0+layerfilter0PL_edgefea_0/bn/cond/switch_t:0 *£
Rlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Qlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Olayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Xlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Zlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Qlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Klayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Tlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Qlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Slayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Qlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Zlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
\layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Slayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Mlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Alayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/decay:0
5layerfilter0PL_edgefea_0/bn/cond/control_dependency:0
*layerfilter0PL_edgefea_0/bn/cond/pred_id:0
+layerfilter0PL_edgefea_0/bn/cond/switch_t:0
glayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
blayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0
ilayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
dlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
-layerfilter0PL_edgefea_0/bn/moments/Squeeze:0
/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1:0X
*layerfilter0PL_edgefea_0/bn/cond/pred_id:0*layerfilter0PL_edgefea_0/bn/cond/pred_id:0Ç
ilayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Zlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1:0\layerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1Ã
glayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0Xlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1¼
dlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0Tlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
-layerfilter0PL_edgefea_0/bn/moments/Squeeze:0Zlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1¸
blayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0Rlayerfilter0PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
ö
,layerfilter0PL_edgefea_0/bn/cond/cond_text_1*layerfilter0PL_edgefea_0/bn/cond/pred_id:0+layerfilter0PL_edgefea_0/bn/cond/switch_f:0*ì
7layerfilter0PL_edgefea_0/bn/cond/control_dependency_1:0
*layerfilter0PL_edgefea_0/bn/cond/pred_id:0
+layerfilter0PL_edgefea_0/bn/cond/switch_f:0X
*layerfilter0PL_edgefea_0/bn/cond/pred_id:0*layerfilter0PL_edgefea_0/bn/cond/pred_id:0
Ë
,layerfilter0PL_edgefea_0/bn/cond_1/cond_text,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0-layerfilter0PL_edgefea_0/bn/cond_1/switch_t:0 *»
4layerfilter0PL_edgefea_0/bn/cond_1/Identity/Switch:1
-layerfilter0PL_edgefea_0/bn/cond_1/Identity:0
6layerfilter0PL_edgefea_0/bn/cond_1/Identity_1/Switch:1
/layerfilter0PL_edgefea_0/bn/cond_1/Identity_1:0
,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0
-layerfilter0PL_edgefea_0/bn/cond_1/switch_t:0
-layerfilter0PL_edgefea_0/bn/moments/Squeeze:0
/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1:0i
/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1:06layerfilter0PL_edgefea_0/bn/cond_1/Identity_1/Switch:1e
-layerfilter0PL_edgefea_0/bn/moments/Squeeze:04layerfilter0PL_edgefea_0/bn/cond_1/Identity/Switch:1\
,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0

.layerfilter0PL_edgefea_0/bn/cond_1/cond_text_1,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0-layerfilter0PL_edgefea_0/bn/cond_1/switch_f:0*
-layerfilter0PL_edgefea_0/bn/cond_1/Switch_1:0
-layerfilter0PL_edgefea_0/bn/cond_1/Switch_1:1
-layerfilter0PL_edgefea_0/bn/cond_1/Switch_2:0
-layerfilter0PL_edgefea_0/bn/cond_1/Switch_2:1
,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0
-layerfilter0PL_edgefea_0/bn/cond_1/switch_f:0
glayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
ilayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
ilayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0-layerfilter0PL_edgefea_0/bn/cond_1/Switch_2:0\
,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0,layerfilter0PL_edgefea_0/bn/cond_1/pred_id:0
glayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0-layerfilter0PL_edgefea_0/bn/cond_1/Switch_1:0
Þ
5layerfilter0PL_self_att_conv_head_0/bn/cond/cond_text5layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:06layerfilter0PL_self_att_conv_head_0/bn/cond/switch_t:0 *³
]layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Zlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
clayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
elayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Vlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
_layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
^layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
\layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
elayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
glayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
^layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Xlayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Llayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
@layerfilter0PL_self_att_conv_head_0/bn/cond/control_dependency:0
5layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:0
6layerfilter0PL_self_att_conv_head_0/bn/cond/switch_t:0
}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
xlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
zlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
8layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1:0Ý
zlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0_layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1¡
8layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze:0elayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1n
5layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:05layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:0è
layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0elayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Ù
xlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0]layerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1¥
:layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1:0glayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1ä
}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0clayerfilter0PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Î
7layerfilter0PL_self_att_conv_head_0/bn/cond/cond_text_15layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:06layerfilter0PL_self_att_conv_head_0/bn/cond/switch_f:0*£
Blayerfilter0PL_self_att_conv_head_0/bn/cond/control_dependency_1:0
5layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:0
6layerfilter0PL_self_att_conv_head_0/bn/cond/switch_f:0n
5layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:05layerfilter0PL_self_att_conv_head_0/bn/cond/pred_id:0

7layerfilter0PL_self_att_conv_head_0/bn/cond_1/cond_text7layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:08layerfilter0PL_self_att_conv_head_0/bn/cond_1/switch_t:0 *Õ
?layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity/Switch:1
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity:0
Alayerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
:layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity_1:0
7layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/switch_t:0
8layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1:0r
7layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:07layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:0{
8layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze:0?layerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity/Switch:1
:layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1:0Alayerfilter0PL_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
ú	
9layerfilter0PL_self_att_conv_head_0/bn/cond_1/cond_text_17layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:08layerfilter0PL_self_att_conv_head_0/bn/cond_1/switch_f:0*É
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_1:0
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_1:1
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_2:0
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_2:1
7layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter0PL_self_att_conv_head_0/bn/cond_1/switch_f:0
}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0»
layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:08layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_2:0r
7layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:07layerfilter0PL_self_att_conv_head_0/bn/cond_1/pred_id:0¹
}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:08layerfilter0PL_self_att_conv_head_0/bn/cond_1/Switch_1:0
Þ
5layerfilter0PL_neib_att_conv_head_0/bn/cond/cond_text5layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:06layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_t:0 *³
]layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Zlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
clayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
elayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Vlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
_layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
^layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
\layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
elayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
glayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
^layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Xlayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Llayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
@layerfilter0PL_neib_att_conv_head_0/bn/cond/control_dependency:0
5layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:0
6layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_t:0
}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
xlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
zlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
8layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0è
layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0elayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¡
8layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze:0elayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ù
xlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0]layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1n
5layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:05layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:0ä
}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0clayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1¥
:layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0glayerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1Ý
zlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0_layerfilter0PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Î
7layerfilter0PL_neib_att_conv_head_0/bn/cond/cond_text_15layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:06layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_f:0*£
Blayerfilter0PL_neib_att_conv_head_0/bn/cond/control_dependency_1:0
5layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:0
6layerfilter0PL_neib_att_conv_head_0/bn/cond/switch_f:0n
5layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:05layerfilter0PL_neib_att_conv_head_0/bn/cond/pred_id:0

7layerfilter0PL_neib_att_conv_head_0/bn/cond_1/cond_text7layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:08layerfilter0PL_neib_att_conv_head_0/bn/cond_1/switch_t:0 *Õ
?layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity:0
Alayerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
:layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity_1:0
7layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/switch_t:0
8layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0
:layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0Alayerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1r
7layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:07layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:0{
8layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze:0?layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1
ú	
9layerfilter0PL_neib_att_conv_head_0/bn/cond_1/cond_text_17layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:08layerfilter0PL_neib_att_conv_head_0/bn/cond_1/switch_f:0*É
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_1:0
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_1:1
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_2:0
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_2:1
7layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter0PL_neib_att_conv_head_0/bn/cond_1/switch_f:0
}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0r
7layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:07layerfilter0PL_neib_att_conv_head_0/bn/cond_1/pred_id:0»
layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:08layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_2:0¹
}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:08layerfilter0PL_neib_att_conv_head_0/bn/cond_1/Switch_1:0
û
gapnet00/bn/cond/cond_textgapnet00/bn/cond/pred_id:0gapnet00/bn/cond/switch_t:0 *¡
Bgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Hgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
;gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Dgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Cgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Lgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Cgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
1gapnet00/bn/cond/ExponentialMovingAverage/decay:0
%gapnet00/bn/cond/control_dependency:0
gapnet00/bn/cond/pred_id:0
gapnet00/bn/cond/switch_t:0
Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Bgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage:0
Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
Dgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage:0
gapnet00/bn/moments/Squeeze:0
gapnet00/bn/moments/Squeeze_1:0k
gapnet00/bn/moments/Squeeze:0Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Dgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage:0Dgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1o
gapnet00/bn/moments/Squeeze_1:0Lgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Bgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage:0Bgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:0Hgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:18
gapnet00/bn/cond/pred_id:0gapnet00/bn/cond/pred_id:0
ö
gapnet00/bn/cond/cond_text_1gapnet00/bn/cond/pred_id:0gapnet00/bn/cond/switch_f:0*
'gapnet00/bn/cond/control_dependency_1:0
gapnet00/bn/cond/pred_id:0
gapnet00/bn/cond/switch_f:08
gapnet00/bn/cond/pred_id:0gapnet00/bn/cond/pred_id:0
»
gapnet00/bn/cond_1/cond_textgapnet00/bn/cond_1/pred_id:0gapnet00/bn/cond_1/switch_t:0 *Û
$gapnet00/bn/cond_1/Identity/Switch:1
gapnet00/bn/cond_1/Identity:0
&gapnet00/bn/cond_1/Identity_1/Switch:1
gapnet00/bn/cond_1/Identity_1:0
gapnet00/bn/cond_1/pred_id:0
gapnet00/bn/cond_1/switch_t:0
gapnet00/bn/moments/Squeeze:0
gapnet00/bn/moments/Squeeze_1:0<
gapnet00/bn/cond_1/pred_id:0gapnet00/bn/cond_1/pred_id:0E
gapnet00/bn/moments/Squeeze:0$gapnet00/bn/cond_1/Identity/Switch:1I
gapnet00/bn/moments/Squeeze_1:0&gapnet00/bn/cond_1/Identity_1/Switch:1
Á
gapnet00/bn/cond_1/cond_text_1gapnet00/bn/cond_1/pred_id:0gapnet00/bn/cond_1/switch_f:0*á
gapnet00/bn/cond_1/Switch_1:0
gapnet00/bn/cond_1/Switch_1:1
gapnet00/bn/cond_1/Switch_2:0
gapnet00/bn/cond_1/Switch_2:1
gapnet00/bn/cond_1/pred_id:0
gapnet00/bn/cond_1/switch_f:0
Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0j
Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet00/bn/cond_1/Switch_2:0<
gapnet00/bn/cond_1/pred_id:0gapnet00/bn/cond_1/pred_id:0h
Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet00/bn/cond_1/Switch_1:0
á
gapnet01PL/bn/cond/cond_textgapnet01PL/bn/cond/pred_id:0gapnet01PL/bn/cond/switch_t:0 *
Dgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Cgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Agapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Jgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Lgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Cgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
=gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Fgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Cgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Egapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Cgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Lgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Ngapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Egapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
?gapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
3gapnet01PL/bn/cond/ExponentialMovingAverage/decay:0
'gapnet01PL/bn/cond/control_dependency:0
gapnet01PL/bn/cond/pred_id:0
gapnet01PL/bn/cond/switch_t:0
Kgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Fgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage:0
Mgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
Hgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage:0
gapnet01PL/bn/moments/Squeeze:0
!gapnet01PL/bn/moments/Squeeze_1:0
Kgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0Jgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Hgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage:0Fgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1s
!gapnet01PL/bn/moments/Squeeze_1:0Ngapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1<
gapnet01PL/bn/cond/pred_id:0gapnet01PL/bn/cond/pred_id:0
Fgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage:0Dgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1o
gapnet01PL/bn/moments/Squeeze:0Lgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Mgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Lgapnet01PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1

gapnet01PL/bn/cond/cond_text_1gapnet01PL/bn/cond/pred_id:0gapnet01PL/bn/cond/switch_f:0*¦
)gapnet01PL/bn/cond/control_dependency_1:0
gapnet01PL/bn/cond/pred_id:0
gapnet01PL/bn/cond/switch_f:0<
gapnet01PL/bn/cond/pred_id:0gapnet01PL/bn/cond/pred_id:0
Ý
gapnet01PL/bn/cond_1/cond_textgapnet01PL/bn/cond_1/pred_id:0gapnet01PL/bn/cond_1/switch_t:0 *÷
&gapnet01PL/bn/cond_1/Identity/Switch:1
gapnet01PL/bn/cond_1/Identity:0
(gapnet01PL/bn/cond_1/Identity_1/Switch:1
!gapnet01PL/bn/cond_1/Identity_1:0
gapnet01PL/bn/cond_1/pred_id:0
gapnet01PL/bn/cond_1/switch_t:0
gapnet01PL/bn/moments/Squeeze:0
!gapnet01PL/bn/moments/Squeeze_1:0@
gapnet01PL/bn/cond_1/pred_id:0gapnet01PL/bn/cond_1/pred_id:0M
!gapnet01PL/bn/moments/Squeeze_1:0(gapnet01PL/bn/cond_1/Identity_1/Switch:1I
gapnet01PL/bn/moments/Squeeze:0&gapnet01PL/bn/cond_1/Identity/Switch:1
ë
 gapnet01PL/bn/cond_1/cond_text_1gapnet01PL/bn/cond_1/pred_id:0gapnet01PL/bn/cond_1/switch_f:0*
gapnet01PL/bn/cond_1/Switch_1:0
gapnet01PL/bn/cond_1/Switch_1:1
gapnet01PL/bn/cond_1/Switch_2:0
gapnet01PL/bn/cond_1/Switch_2:1
gapnet01PL/bn/cond_1/pred_id:0
gapnet01PL/bn/cond_1/switch_f:0
Kgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Mgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0@
gapnet01PL/bn/cond_1/pred_id:0gapnet01PL/bn/cond_1/pred_id:0n
Kgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet01PL/bn/cond_1/Switch_1:0p
Mgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet01PL/bn/cond_1/Switch_2:0
ø
3layerfilter1PL_newfea_conv_head_0/bn/cond/cond_text3layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:04layerfilter1PL_newfea_conv_head_0/bn/cond/switch_t:0 *Ó
[layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter1PL_newfea_conv_head_0/bn/cond/control_dependency:0
3layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:0
4layerfilter1PL_newfea_conv_head_0/bn/cond/switch_t:0
ylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze:0
8layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1:0×
vlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1¡
8layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1:0elayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1Ó
tlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
6layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze:0clayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1â
{layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1j
3layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:03layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:0Þ
ylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter1PL_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
¾
5layerfilter1PL_newfea_conv_head_0/bn/cond/cond_text_13layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:04layerfilter1PL_newfea_conv_head_0/bn/cond/switch_f:0*
@layerfilter1PL_newfea_conv_head_0/bn/cond/control_dependency_1:0
3layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:0
4layerfilter1PL_newfea_conv_head_0/bn/cond/switch_f:0j
3layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:03layerfilter1PL_newfea_conv_head_0/bn/cond/pred_id:0
ä
5layerfilter1PL_newfea_conv_head_0/bn/cond_1/cond_text5layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:06layerfilter1PL_newfea_conv_head_0/bn/cond_1/switch_t:0 *¹
=layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity/Switch:1
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity:0
?layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1
8layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity_1:0
5layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:0
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/switch_t:0
6layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze:0
8layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1:0{
8layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1:0?layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1w
6layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze:0=layerfilter1PL_newfea_conv_head_0/bn/cond_1/Identity/Switch:1n
5layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:05layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:0
Ð	
7layerfilter1PL_newfea_conv_head_0/bn/cond_1/cond_text_15layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:06layerfilter1PL_newfea_conv_head_0/bn/cond_1/switch_f:0*¥
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_1:0
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_1:1
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_2:0
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_2:1
5layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:0
6layerfilter1PL_newfea_conv_head_0/bn/cond_1/switch_f:0
ylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_1:0n
5layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:05layerfilter1PL_newfea_conv_head_0/bn/cond_1/pred_id:0µ
{layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter1PL_newfea_conv_head_0/bn/cond_1/Switch_2:0
­
*layerfilter1PL_edgefea_0/bn/cond/cond_text*layerfilter1PL_edgefea_0/bn/cond/pred_id:0+layerfilter1PL_edgefea_0/bn/cond/switch_t:0 *£
Rlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Qlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Olayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Xlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Zlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Qlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Klayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Tlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Qlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Slayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Qlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Zlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
\layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Slayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Mlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Alayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/decay:0
5layerfilter1PL_edgefea_0/bn/cond/control_dependency:0
*layerfilter1PL_edgefea_0/bn/cond/pred_id:0
+layerfilter1PL_edgefea_0/bn/cond/switch_t:0
glayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
blayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0
ilayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
dlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
-layerfilter1PL_edgefea_0/bn/moments/Squeeze:0
/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1:0
-layerfilter1PL_edgefea_0/bn/moments/Squeeze:0Zlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ã
glayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0Xlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1:0\layerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1¼
dlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0Tlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1X
*layerfilter1PL_edgefea_0/bn/cond/pred_id:0*layerfilter1PL_edgefea_0/bn/cond/pred_id:0Ç
ilayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Zlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¸
blayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0Rlayerfilter1PL_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
ö
,layerfilter1PL_edgefea_0/bn/cond/cond_text_1*layerfilter1PL_edgefea_0/bn/cond/pred_id:0+layerfilter1PL_edgefea_0/bn/cond/switch_f:0*ì
7layerfilter1PL_edgefea_0/bn/cond/control_dependency_1:0
*layerfilter1PL_edgefea_0/bn/cond/pred_id:0
+layerfilter1PL_edgefea_0/bn/cond/switch_f:0X
*layerfilter1PL_edgefea_0/bn/cond/pred_id:0*layerfilter1PL_edgefea_0/bn/cond/pred_id:0
Ë
,layerfilter1PL_edgefea_0/bn/cond_1/cond_text,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0-layerfilter1PL_edgefea_0/bn/cond_1/switch_t:0 *»
4layerfilter1PL_edgefea_0/bn/cond_1/Identity/Switch:1
-layerfilter1PL_edgefea_0/bn/cond_1/Identity:0
6layerfilter1PL_edgefea_0/bn/cond_1/Identity_1/Switch:1
/layerfilter1PL_edgefea_0/bn/cond_1/Identity_1:0
,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0
-layerfilter1PL_edgefea_0/bn/cond_1/switch_t:0
-layerfilter1PL_edgefea_0/bn/moments/Squeeze:0
/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1:0e
-layerfilter1PL_edgefea_0/bn/moments/Squeeze:04layerfilter1PL_edgefea_0/bn/cond_1/Identity/Switch:1i
/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1:06layerfilter1PL_edgefea_0/bn/cond_1/Identity_1/Switch:1\
,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0

.layerfilter1PL_edgefea_0/bn/cond_1/cond_text_1,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0-layerfilter1PL_edgefea_0/bn/cond_1/switch_f:0*
-layerfilter1PL_edgefea_0/bn/cond_1/Switch_1:0
-layerfilter1PL_edgefea_0/bn/cond_1/Switch_1:1
-layerfilter1PL_edgefea_0/bn/cond_1/Switch_2:0
-layerfilter1PL_edgefea_0/bn/cond_1/Switch_2:1
,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0
-layerfilter1PL_edgefea_0/bn/cond_1/switch_f:0
glayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
ilayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
ilayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0-layerfilter1PL_edgefea_0/bn/cond_1/Switch_2:0\
,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0,layerfilter1PL_edgefea_0/bn/cond_1/pred_id:0
glayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0-layerfilter1PL_edgefea_0/bn/cond_1/Switch_1:0
Þ
5layerfilter1PL_self_att_conv_head_0/bn/cond/cond_text5layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:06layerfilter1PL_self_att_conv_head_0/bn/cond/switch_t:0 *³
]layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Zlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
clayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
elayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Vlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
_layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
^layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
\layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
elayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
glayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
^layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Xlayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Llayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
@layerfilter1PL_self_att_conv_head_0/bn/cond/control_dependency:0
5layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:0
6layerfilter1PL_self_att_conv_head_0/bn/cond/switch_t:0
}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
xlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
zlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
8layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1:0¥
:layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1:0glayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1¡
8layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze:0elayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1è
layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0elayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1ä
}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0clayerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Ý
zlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0_layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1n
5layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:05layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:0Ù
xlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0]layerfilter1PL_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Î
7layerfilter1PL_self_att_conv_head_0/bn/cond/cond_text_15layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:06layerfilter1PL_self_att_conv_head_0/bn/cond/switch_f:0*£
Blayerfilter1PL_self_att_conv_head_0/bn/cond/control_dependency_1:0
5layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:0
6layerfilter1PL_self_att_conv_head_0/bn/cond/switch_f:0n
5layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:05layerfilter1PL_self_att_conv_head_0/bn/cond/pred_id:0

7layerfilter1PL_self_att_conv_head_0/bn/cond_1/cond_text7layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:08layerfilter1PL_self_att_conv_head_0/bn/cond_1/switch_t:0 *Õ
?layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity/Switch:1
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity:0
Alayerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
:layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity_1:0
7layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/switch_t:0
8layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1:0{
8layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze:0?layerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity/Switch:1r
7layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:07layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:0
:layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1:0Alayerfilter1PL_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
ú	
9layerfilter1PL_self_att_conv_head_0/bn/cond_1/cond_text_17layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:08layerfilter1PL_self_att_conv_head_0/bn/cond_1/switch_f:0*É
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_1:0
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_1:1
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_2:0
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_2:1
7layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter1PL_self_att_conv_head_0/bn/cond_1/switch_f:0
}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0¹
}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:08layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_1:0»
layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:08layerfilter1PL_self_att_conv_head_0/bn/cond_1/Switch_2:0r
7layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:07layerfilter1PL_self_att_conv_head_0/bn/cond_1/pred_id:0
Þ
5layerfilter1PL_neib_att_conv_head_0/bn/cond/cond_text5layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:06layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_t:0 *³
]layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Zlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
clayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
elayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Vlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
_layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
^layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
\layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
elayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
glayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
^layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Xlayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Llayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
@layerfilter1PL_neib_att_conv_head_0/bn/cond/control_dependency:0
5layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:0
6layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_t:0
}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
xlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
zlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
8layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0¡
8layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze:0elayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1ä
}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0clayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Ý
zlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0_layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Ù
xlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0]layerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1n
5layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:05layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:0¥
:layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0glayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1è
layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0elayerfilter1PL_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Î
7layerfilter1PL_neib_att_conv_head_0/bn/cond/cond_text_15layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:06layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_f:0*£
Blayerfilter1PL_neib_att_conv_head_0/bn/cond/control_dependency_1:0
5layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:0
6layerfilter1PL_neib_att_conv_head_0/bn/cond/switch_f:0n
5layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:05layerfilter1PL_neib_att_conv_head_0/bn/cond/pred_id:0

7layerfilter1PL_neib_att_conv_head_0/bn/cond_1/cond_text7layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:08layerfilter1PL_neib_att_conv_head_0/bn/cond_1/switch_t:0 *Õ
?layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity:0
Alayerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
:layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity_1:0
7layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/switch_t:0
8layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze:0
:layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0{
8layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze:0?layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1r
7layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:07layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:0
:layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1:0Alayerfilter1PL_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
ú	
9layerfilter1PL_neib_att_conv_head_0/bn/cond_1/cond_text_17layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:08layerfilter1PL_neib_att_conv_head_0/bn/cond_1/switch_f:0*É
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_1:0
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_1:1
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_2:0
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_2:1
7layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:0
8layerfilter1PL_neib_att_conv_head_0/bn/cond_1/switch_f:0
}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0r
7layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:07layerfilter1PL_neib_att_conv_head_0/bn/cond_1/pred_id:0»
layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:08layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_2:0¹
}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:08layerfilter1PL_neib_att_conv_head_0/bn/cond_1/Switch_1:0
û
gapnet10/bn/cond/cond_textgapnet10/bn/cond/pred_id:0gapnet10/bn/cond/switch_t:0 *¡
Bgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Hgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
;gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Dgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Cgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Lgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Cgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
1gapnet10/bn/cond/ExponentialMovingAverage/decay:0
%gapnet10/bn/cond/control_dependency:0
gapnet10/bn/cond/pred_id:0
gapnet10/bn/cond/switch_t:0
Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Bgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage:0
Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
Dgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage:0
gapnet10/bn/moments/Squeeze:0
gapnet10/bn/moments/Squeeze_1:0
Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:18
gapnet10/bn/cond/pred_id:0gapnet10/bn/cond/pred_id:0k
gapnet10/bn/moments/Squeeze:0Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:0Hgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Dgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage:0Dgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Bgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage:0Bgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1o
gapnet10/bn/moments/Squeeze_1:0Lgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
ö
gapnet10/bn/cond/cond_text_1gapnet10/bn/cond/pred_id:0gapnet10/bn/cond/switch_f:0*
'gapnet10/bn/cond/control_dependency_1:0
gapnet10/bn/cond/pred_id:0
gapnet10/bn/cond/switch_f:08
gapnet10/bn/cond/pred_id:0gapnet10/bn/cond/pred_id:0
»
gapnet10/bn/cond_1/cond_textgapnet10/bn/cond_1/pred_id:0gapnet10/bn/cond_1/switch_t:0 *Û
$gapnet10/bn/cond_1/Identity/Switch:1
gapnet10/bn/cond_1/Identity:0
&gapnet10/bn/cond_1/Identity_1/Switch:1
gapnet10/bn/cond_1/Identity_1:0
gapnet10/bn/cond_1/pred_id:0
gapnet10/bn/cond_1/switch_t:0
gapnet10/bn/moments/Squeeze:0
gapnet10/bn/moments/Squeeze_1:0E
gapnet10/bn/moments/Squeeze:0$gapnet10/bn/cond_1/Identity/Switch:1<
gapnet10/bn/cond_1/pred_id:0gapnet10/bn/cond_1/pred_id:0I
gapnet10/bn/moments/Squeeze_1:0&gapnet10/bn/cond_1/Identity_1/Switch:1
Á
gapnet10/bn/cond_1/cond_text_1gapnet10/bn/cond_1/pred_id:0gapnet10/bn/cond_1/switch_f:0*á
gapnet10/bn/cond_1/Switch_1:0
gapnet10/bn/cond_1/Switch_1:1
gapnet10/bn/cond_1/Switch_2:0
gapnet10/bn/cond_1/Switch_2:1
gapnet10/bn/cond_1/pred_id:0
gapnet10/bn/cond_1/switch_f:0
Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0<
gapnet10/bn/cond_1/pred_id:0gapnet10/bn/cond_1/pred_id:0h
Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet10/bn/cond_1/Switch_1:0j
Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet10/bn/cond_1/Switch_2:0
á
gapnet11PL/bn/cond/cond_textgapnet11PL/bn/cond/pred_id:0gapnet11PL/bn/cond/switch_t:0 *
Dgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Cgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Agapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Jgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Lgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Cgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
=gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Fgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Cgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Egapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Cgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Lgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Ngapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Egapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
?gapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
3gapnet11PL/bn/cond/ExponentialMovingAverage/decay:0
'gapnet11PL/bn/cond/control_dependency:0
gapnet11PL/bn/cond/pred_id:0
gapnet11PL/bn/cond/switch_t:0
Kgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Fgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage:0
Mgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
Hgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage:0
gapnet11PL/bn/moments/Squeeze:0
!gapnet11PL/bn/moments/Squeeze_1:0
Hgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage:0Fgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Mgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Lgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1o
gapnet11PL/bn/moments/Squeeze:0Lgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Fgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage:0Dgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1<
gapnet11PL/bn/cond/pred_id:0gapnet11PL/bn/cond/pred_id:0
Kgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0Jgapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1s
!gapnet11PL/bn/moments/Squeeze_1:0Ngapnet11PL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1

gapnet11PL/bn/cond/cond_text_1gapnet11PL/bn/cond/pred_id:0gapnet11PL/bn/cond/switch_f:0*¦
)gapnet11PL/bn/cond/control_dependency_1:0
gapnet11PL/bn/cond/pred_id:0
gapnet11PL/bn/cond/switch_f:0<
gapnet11PL/bn/cond/pred_id:0gapnet11PL/bn/cond/pred_id:0
Ý
gapnet11PL/bn/cond_1/cond_textgapnet11PL/bn/cond_1/pred_id:0gapnet11PL/bn/cond_1/switch_t:0 *÷
&gapnet11PL/bn/cond_1/Identity/Switch:1
gapnet11PL/bn/cond_1/Identity:0
(gapnet11PL/bn/cond_1/Identity_1/Switch:1
!gapnet11PL/bn/cond_1/Identity_1:0
gapnet11PL/bn/cond_1/pred_id:0
gapnet11PL/bn/cond_1/switch_t:0
gapnet11PL/bn/moments/Squeeze:0
!gapnet11PL/bn/moments/Squeeze_1:0M
!gapnet11PL/bn/moments/Squeeze_1:0(gapnet11PL/bn/cond_1/Identity_1/Switch:1@
gapnet11PL/bn/cond_1/pred_id:0gapnet11PL/bn/cond_1/pred_id:0I
gapnet11PL/bn/moments/Squeeze:0&gapnet11PL/bn/cond_1/Identity/Switch:1
ë
 gapnet11PL/bn/cond_1/cond_text_1gapnet11PL/bn/cond_1/pred_id:0gapnet11PL/bn/cond_1/switch_f:0*
gapnet11PL/bn/cond_1/Switch_1:0
gapnet11PL/bn/cond_1/Switch_1:1
gapnet11PL/bn/cond_1/Switch_2:0
gapnet11PL/bn/cond_1/Switch_2:1
gapnet11PL/bn/cond_1/pred_id:0
gapnet11PL/bn/cond_1/switch_f:0
Kgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Mgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0n
Kgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet11PL/bn/cond_1/Switch_1:0@
gapnet11PL/bn/cond_1/pred_id:0gapnet11PL/bn/cond_1/pred_id:0p
Mgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet11PL/bn/cond_1/Switch_2:0
á
aggPL/bn/cond/cond_textaggPL/bn/cond/pred_id:0aggPL/bn/cond/switch_t:0 *
AaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/read:0
<aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage:0
CaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
>aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage:0
?aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
<aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
EaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
GaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
8aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
AaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
@aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
>aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
GaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
IaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
@aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
:aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
.aggPL/bn/cond/ExponentialMovingAverage/decay:0
"aggPL/bn/cond/control_dependency:0
aggPL/bn/cond/pred_id:0
aggPL/bn/cond/switch_t:0
aggPL/bn/moments/Squeeze:0
aggPL/bn/moments/Squeeze_1:0
CaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0GaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
<aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage:0?aggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:12
aggPL/bn/cond/pred_id:0aggPL/bn/cond/pred_id:0
AaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/read:0EaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1e
aggPL/bn/moments/Squeeze:0GaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1i
aggPL/bn/moments/Squeeze_1:0IaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
>aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage:0AaggPL/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Þ
aggPL/bn/cond/cond_text_1aggPL/bn/cond/pred_id:0aggPL/bn/cond/switch_f:0*
$aggPL/bn/cond/control_dependency_1:0
aggPL/bn/cond/pred_id:0
aggPL/bn/cond/switch_f:02
aggPL/bn/cond/pred_id:0aggPL/bn/cond/pred_id:0

aggPL/bn/cond_1/cond_textaggPL/bn/cond_1/pred_id:0aggPL/bn/cond_1/switch_t:0 *±
!aggPL/bn/cond_1/Identity/Switch:1
aggPL/bn/cond_1/Identity:0
#aggPL/bn/cond_1/Identity_1/Switch:1
aggPL/bn/cond_1/Identity_1:0
aggPL/bn/cond_1/pred_id:0
aggPL/bn/cond_1/switch_t:0
aggPL/bn/moments/Squeeze:0
aggPL/bn/moments/Squeeze_1:0C
aggPL/bn/moments/Squeeze_1:0#aggPL/bn/cond_1/Identity_1/Switch:16
aggPL/bn/cond_1/pred_id:0aggPL/bn/cond_1/pred_id:0?
aggPL/bn/moments/Squeeze:0!aggPL/bn/cond_1/Identity/Switch:1

aggPL/bn/cond_1/cond_text_1aggPL/bn/cond_1/pred_id:0aggPL/bn/cond_1/switch_f:0*«
AaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/read:0
CaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
aggPL/bn/cond_1/Switch_1:0
aggPL/bn/cond_1/Switch_1:1
aggPL/bn/cond_1/Switch_2:0
aggPL/bn/cond_1/Switch_2:1
aggPL/bn/cond_1/pred_id:0
aggPL/bn/cond_1/switch_f:0_
AaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/read:0aggPL/bn/cond_1/Switch_1:06
aggPL/bn/cond_1/pred_id:0aggPL/bn/cond_1/pred_id:0a
CaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0aggPL/bn/cond_1/Switch_2:0
á
PLfc1/bn/cond/cond_textPLfc1/bn/cond/pred_id:0PLfc1/bn/cond/switch_t:0 *
APLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
<PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage:0
CPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
>PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage:0
?PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
<PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
EPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
GPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
8PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
APLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
@PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
>PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
GPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
IPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
@PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
:PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
.PLfc1/bn/cond/ExponentialMovingAverage/decay:0
"PLfc1/bn/cond/control_dependency:0
PLfc1/bn/cond/pred_id:0
PLfc1/bn/cond/switch_t:0
PLfc1/bn/moments/Squeeze:0
PLfc1/bn/moments/Squeeze_1:0
>PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage:0APLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
<PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage:0?PLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:12
PLfc1/bn/cond/pred_id:0PLfc1/bn/cond/pred_id:0i
PLfc1/bn/moments/Squeeze_1:0IPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
CPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0GPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1e
PLfc1/bn/moments/Squeeze:0GPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
APLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/read:0EPLfc1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Þ
PLfc1/bn/cond/cond_text_1PLfc1/bn/cond/pred_id:0PLfc1/bn/cond/switch_f:0*
$PLfc1/bn/cond/control_dependency_1:0
PLfc1/bn/cond/pred_id:0
PLfc1/bn/cond/switch_f:02
PLfc1/bn/cond/pred_id:0PLfc1/bn/cond/pred_id:0

PLfc1/bn/cond_1/cond_textPLfc1/bn/cond_1/pred_id:0PLfc1/bn/cond_1/switch_t:0 *±
!PLfc1/bn/cond_1/Identity/Switch:1
PLfc1/bn/cond_1/Identity:0
#PLfc1/bn/cond_1/Identity_1/Switch:1
PLfc1/bn/cond_1/Identity_1:0
PLfc1/bn/cond_1/pred_id:0
PLfc1/bn/cond_1/switch_t:0
PLfc1/bn/moments/Squeeze:0
PLfc1/bn/moments/Squeeze_1:0?
PLfc1/bn/moments/Squeeze:0!PLfc1/bn/cond_1/Identity/Switch:16
PLfc1/bn/cond_1/pred_id:0PLfc1/bn/cond_1/pred_id:0C
PLfc1/bn/moments/Squeeze_1:0#PLfc1/bn/cond_1/Identity_1/Switch:1

PLfc1/bn/cond_1/cond_text_1PLfc1/bn/cond_1/pred_id:0PLfc1/bn/cond_1/switch_f:0*«
APLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
CPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
PLfc1/bn/cond_1/Switch_1:0
PLfc1/bn/cond_1/Switch_1:1
PLfc1/bn/cond_1/Switch_2:0
PLfc1/bn/cond_1/Switch_2:1
PLfc1/bn/cond_1/pred_id:0
PLfc1/bn/cond_1/switch_f:0a
CPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0PLfc1/bn/cond_1/Switch_2:0_
APLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/read:0PLfc1/bn/cond_1/Switch_1:06
PLfc1/bn/cond_1/pred_id:0PLfc1/bn/cond_1/pred_id:0

PLdp1/cond/cond_textPLdp1/cond/pred_id:0PLdp1/cond/switch_t:0 *Ö
PLdp1/cond/dropout/Cast:0
!PLdp1/cond/dropout/GreaterEqual:0
PLdp1/cond/dropout/Shape:0
PLdp1/cond/dropout/mul/Switch:1
PLdp1/cond/dropout/mul:0
PLdp1/cond/dropout/mul_1:0
1PLdp1/cond/dropout/random_uniform/RandomUniform:0
'PLdp1/cond/dropout/random_uniform/max:0
'PLdp1/cond/dropout/random_uniform/min:0
'PLdp1/cond/dropout/random_uniform/mul:0
'PLdp1/cond/dropout/random_uniform/sub:0
#PLdp1/cond/dropout/random_uniform:0
PLdp1/cond/dropout/rate:0
PLdp1/cond/dropout/sub/x:0
PLdp1/cond/dropout/sub:0
PLdp1/cond/dropout/truediv/x:0
PLdp1/cond/dropout/truediv:0
PLdp1/cond/pred_id:0
PLdp1/cond/switch_t:0
PLfc1/Relu:0/
PLfc1/Relu:0PLdp1/cond/dropout/mul/Switch:1,
PLdp1/cond/pred_id:0PLdp1/cond/pred_id:0

PLdp1/cond/cond_text_1PLdp1/cond/pred_id:0PLdp1/cond/switch_f:0*¾
PLdp1/cond/Switch_1:0
PLdp1/cond/Switch_1:1
PLdp1/cond/pred_id:0
PLdp1/cond/switch_f:0
PLfc1/Relu:0%
PLfc1/Relu:0PLdp1/cond/Switch_1:0,
PLdp1/cond/pred_id:0PLdp1/cond/pred_id:0
á
PLfc2/bn/cond/cond_textPLfc2/bn/cond/pred_id:0PLfc2/bn/cond/switch_t:0 *
APLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/read:0
<PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage:0
CPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
>PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage:0
?PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
<PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
EPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
GPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
8PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
APLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
@PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
>PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
GPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
IPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
@PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
:PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
.PLfc2/bn/cond/ExponentialMovingAverage/decay:0
"PLfc2/bn/cond/control_dependency:0
PLfc2/bn/cond/pred_id:0
PLfc2/bn/cond/switch_t:0
PLfc2/bn/moments/Squeeze:0
PLfc2/bn/moments/Squeeze_1:0e
PLfc2/bn/moments/Squeeze:0GPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
>PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage:0APLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1i
PLfc2/bn/moments/Squeeze_1:0IPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
CPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0GPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
<PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage:0?PLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
APLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/read:0EPLfc2/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:12
PLfc2/bn/cond/pred_id:0PLfc2/bn/cond/pred_id:0
Þ
PLfc2/bn/cond/cond_text_1PLfc2/bn/cond/pred_id:0PLfc2/bn/cond/switch_f:0*
$PLfc2/bn/cond/control_dependency_1:0
PLfc2/bn/cond/pred_id:0
PLfc2/bn/cond/switch_f:02
PLfc2/bn/cond/pred_id:0PLfc2/bn/cond/pred_id:0

PLfc2/bn/cond_1/cond_textPLfc2/bn/cond_1/pred_id:0PLfc2/bn/cond_1/switch_t:0 *±
!PLfc2/bn/cond_1/Identity/Switch:1
PLfc2/bn/cond_1/Identity:0
#PLfc2/bn/cond_1/Identity_1/Switch:1
PLfc2/bn/cond_1/Identity_1:0
PLfc2/bn/cond_1/pred_id:0
PLfc2/bn/cond_1/switch_t:0
PLfc2/bn/moments/Squeeze:0
PLfc2/bn/moments/Squeeze_1:0C
PLfc2/bn/moments/Squeeze_1:0#PLfc2/bn/cond_1/Identity_1/Switch:16
PLfc2/bn/cond_1/pred_id:0PLfc2/bn/cond_1/pred_id:0?
PLfc2/bn/moments/Squeeze:0!PLfc2/bn/cond_1/Identity/Switch:1

PLfc2/bn/cond_1/cond_text_1PLfc2/bn/cond_1/pred_id:0PLfc2/bn/cond_1/switch_f:0*«
APLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/read:0
CPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
PLfc2/bn/cond_1/Switch_1:0
PLfc2/bn/cond_1/Switch_1:1
PLfc2/bn/cond_1/Switch_2:0
PLfc2/bn/cond_1/Switch_2:1
PLfc2/bn/cond_1/pred_id:0
PLfc2/bn/cond_1/switch_f:0a
CPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0PLfc2/bn/cond_1/Switch_2:06
PLfc2/bn/cond_1/pred_id:0PLfc2/bn/cond_1/pred_id:0_
APLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/read:0PLfc2/bn/cond_1/Switch_1:0

PLdp2/cond/cond_textPLdp2/cond/pred_id:0PLdp2/cond/switch_t:0 *Ö
PLdp2/cond/dropout/Cast:0
!PLdp2/cond/dropout/GreaterEqual:0
PLdp2/cond/dropout/Shape:0
PLdp2/cond/dropout/mul/Switch:1
PLdp2/cond/dropout/mul:0
PLdp2/cond/dropout/mul_1:0
1PLdp2/cond/dropout/random_uniform/RandomUniform:0
'PLdp2/cond/dropout/random_uniform/max:0
'PLdp2/cond/dropout/random_uniform/min:0
'PLdp2/cond/dropout/random_uniform/mul:0
'PLdp2/cond/dropout/random_uniform/sub:0
#PLdp2/cond/dropout/random_uniform:0
PLdp2/cond/dropout/rate:0
PLdp2/cond/dropout/sub/x:0
PLdp2/cond/dropout/sub:0
PLdp2/cond/dropout/truediv/x:0
PLdp2/cond/dropout/truediv:0
PLdp2/cond/pred_id:0
PLdp2/cond/switch_t:0
PLfc2/Relu:0/
PLfc2/Relu:0PLdp2/cond/dropout/mul/Switch:1,
PLdp2/cond/pred_id:0PLdp2/cond/pred_id:0

PLdp2/cond/cond_text_1PLdp2/cond/pred_id:0PLdp2/cond/switch_f:0*¾
PLdp2/cond/Switch_1:0
PLdp2/cond/Switch_1:1
PLdp2/cond/pred_id:0
PLdp2/cond/switch_f:0
PLfc2/Relu:0%
PLfc2/Relu:0PLdp2/cond/Switch_1:0,
PLdp2/cond/pred_id:0PLdp2/cond/pred_id:0"¦©
	variables©©
Û
+layerfilter0PL_newfea_conv_head_0/weights:00layerfilter0PL_newfea_conv_head_0/weights/Assign0layerfilter0PL_newfea_conv_head_0/weights/read:02Flayerfilter0PL_newfea_conv_head_0/weights/Initializer/random_uniform:08
Á
+layerfilter0PL_newfea_conv_head_0/bn/beta:00layerfilter0PL_newfea_conv_head_0/bn/beta/Assign0layerfilter0PL_newfea_conv_head_0/bn/beta/read:02,layerfilter0PL_newfea_conv_head_0/bn/Const:08
Æ
,layerfilter0PL_newfea_conv_head_0/bn/gamma:01layerfilter0PL_newfea_conv_head_0/bn/gamma/Assign1layerfilter0PL_newfea_conv_head_0/bn/gamma/read:02.layerfilter0PL_newfea_conv_head_0/bn/Const_1:08
õ
tlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter0PL_newfea_conv_head_0/bn/layerfilter0PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
·
"layerfilter0PL_edgefea_0/weights:0'layerfilter0PL_edgefea_0/weights/Assign'layerfilter0PL_edgefea_0/weights/read:02=layerfilter0PL_edgefea_0/weights/Initializer/random_uniform:08
ª
!layerfilter0PL_edgefea_0/biases:0&layerfilter0PL_edgefea_0/biases/Assign&layerfilter0PL_edgefea_0/biases/read:023layerfilter0PL_edgefea_0/biases/Initializer/Const:08

"layerfilter0PL_edgefea_0/bn/beta:0'layerfilter0PL_edgefea_0/bn/beta/Assign'layerfilter0PL_edgefea_0/bn/beta/read:02#layerfilter0PL_edgefea_0/bn/Const:08
¢
#layerfilter0PL_edgefea_0/bn/gamma:0(layerfilter0PL_edgefea_0/bn/gamma/Assign(layerfilter0PL_edgefea_0/bn/gamma/read:02%layerfilter0PL_edgefea_0/bn/Const_1:08
¬
blayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0glayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignglayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02tlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
dlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0ilayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignilayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02vlayerfilter0PL_edgefea_0/bn/layerfilter0PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
ã
-layerfilter0PL_self_att_conv_head_0/weights:02layerfilter0PL_self_att_conv_head_0/weights/Assign2layerfilter0PL_self_att_conv_head_0/weights/read:02Hlayerfilter0PL_self_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter0PL_self_att_conv_head_0/biases:01layerfilter0PL_self_att_conv_head_0/biases/Assign1layerfilter0PL_self_att_conv_head_0/biases/read:02>layerfilter0PL_self_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter0PL_self_att_conv_head_0/bn/beta:02layerfilter0PL_self_att_conv_head_0/bn/beta/Assign2layerfilter0PL_self_att_conv_head_0/bn/beta/read:02.layerfilter0PL_self_att_conv_head_0/bn/Const:08
Î
.layerfilter0PL_self_att_conv_head_0/bn/gamma:03layerfilter0PL_self_att_conv_head_0/bn/gamma/Assign3layerfilter0PL_self_att_conv_head_0/bn/gamma/read:020layerfilter0PL_self_att_conv_head_0/bn/Const_1:08

xlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assign}layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

zlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignlayerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter0PL_self_att_conv_head_0/bn/layerfilter0PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
ã
-layerfilter0PL_neib_att_conv_head_0/weights:02layerfilter0PL_neib_att_conv_head_0/weights/Assign2layerfilter0PL_neib_att_conv_head_0/weights/read:02Hlayerfilter0PL_neib_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter0PL_neib_att_conv_head_0/biases:01layerfilter0PL_neib_att_conv_head_0/biases/Assign1layerfilter0PL_neib_att_conv_head_0/biases/read:02>layerfilter0PL_neib_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter0PL_neib_att_conv_head_0/bn/beta:02layerfilter0PL_neib_att_conv_head_0/bn/beta/Assign2layerfilter0PL_neib_att_conv_head_0/bn/beta/read:02.layerfilter0PL_neib_att_conv_head_0/bn/Const:08
Î
.layerfilter0PL_neib_att_conv_head_0/bn/gamma:03layerfilter0PL_neib_att_conv_head_0/bn/gamma/Assign3layerfilter0PL_neib_att_conv_head_0/bn/gamma/read:020layerfilter0PL_neib_att_conv_head_0/bn/Const_1:08

xlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assign}layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

zlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignlayerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter0PL_neib_att_conv_head_0/bn/layerfilter0PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
¢
layerfilter0PL/BiasAdd/biases:0$layerfilter0PL/BiasAdd/biases/Assign$layerfilter0PL/BiasAdd/biases/read:021layerfilter0PL/BiasAdd/biases/Initializer/zeros:08
w
gapnet00/weights:0gapnet00/weights/Assigngapnet00/weights/read:02-gapnet00/weights/Initializer/random_uniform:08
j
gapnet00/biases:0gapnet00/biases/Assigngapnet00/biases/read:02#gapnet00/biases/Initializer/Const:08
]
gapnet00/bn/beta:0gapnet00/bn/beta/Assigngapnet00/bn/beta/read:02gapnet00/bn/Const:08
b
gapnet00/bn/gamma:0gapnet00/bn/gamma/Assigngapnet00/bn/gamma/read:02gapnet00/bn/Const_1:08
¬
Bgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage:0Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/AssignGgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:02Tgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
Dgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage:0Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignIgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Vgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0

gapnet01PL/weights:0gapnet01PL/weights/Assigngapnet01PL/weights/read:02/gapnet01PL/weights/Initializer/random_uniform:08
r
gapnet01PL/biases:0gapnet01PL/biases/Assigngapnet01PL/biases/read:02%gapnet01PL/biases/Initializer/Const:08
e
gapnet01PL/bn/beta:0gapnet01PL/bn/beta/Assigngapnet01PL/bn/beta/read:02gapnet01PL/bn/Const:08
j
gapnet01PL/bn/gamma:0gapnet01PL/bn/gamma/Assigngapnet01PL/bn/gamma/read:02gapnet01PL/bn/Const_1:08
¼
Fgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage:0Kgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/AssignKgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/read:02Xgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
Ä
Hgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage:0Mgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignMgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Zgapnet01PL/bn/gapnet01PL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Û
+layerfilter1PL_newfea_conv_head_0/weights:00layerfilter1PL_newfea_conv_head_0/weights/Assign0layerfilter1PL_newfea_conv_head_0/weights/read:02Flayerfilter1PL_newfea_conv_head_0/weights/Initializer/random_uniform:08
Á
+layerfilter1PL_newfea_conv_head_0/bn/beta:00layerfilter1PL_newfea_conv_head_0/bn/beta/Assign0layerfilter1PL_newfea_conv_head_0/bn/beta/read:02,layerfilter1PL_newfea_conv_head_0/bn/Const:08
Æ
,layerfilter1PL_newfea_conv_head_0/bn/gamma:01layerfilter1PL_newfea_conv_head_0/bn/gamma/Assign1layerfilter1PL_newfea_conv_head_0/bn/gamma/read:02.layerfilter1PL_newfea_conv_head_0/bn/Const_1:08
õ
tlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1PL_newfea_conv_head_0/bn/layerfilter1PL_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
·
"layerfilter1PL_edgefea_0/weights:0'layerfilter1PL_edgefea_0/weights/Assign'layerfilter1PL_edgefea_0/weights/read:02=layerfilter1PL_edgefea_0/weights/Initializer/random_uniform:08
ª
!layerfilter1PL_edgefea_0/biases:0&layerfilter1PL_edgefea_0/biases/Assign&layerfilter1PL_edgefea_0/biases/read:023layerfilter1PL_edgefea_0/biases/Initializer/Const:08

"layerfilter1PL_edgefea_0/bn/beta:0'layerfilter1PL_edgefea_0/bn/beta/Assign'layerfilter1PL_edgefea_0/bn/beta/read:02#layerfilter1PL_edgefea_0/bn/Const:08
¢
#layerfilter1PL_edgefea_0/bn/gamma:0(layerfilter1PL_edgefea_0/bn/gamma/Assign(layerfilter1PL_edgefea_0/bn/gamma/read:02%layerfilter1PL_edgefea_0/bn/Const_1:08
¬
blayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0glayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignglayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02tlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
dlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0ilayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignilayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02vlayerfilter1PL_edgefea_0/bn/layerfilter1PL_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
ã
-layerfilter1PL_self_att_conv_head_0/weights:02layerfilter1PL_self_att_conv_head_0/weights/Assign2layerfilter1PL_self_att_conv_head_0/weights/read:02Hlayerfilter1PL_self_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter1PL_self_att_conv_head_0/biases:01layerfilter1PL_self_att_conv_head_0/biases/Assign1layerfilter1PL_self_att_conv_head_0/biases/read:02>layerfilter1PL_self_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter1PL_self_att_conv_head_0/bn/beta:02layerfilter1PL_self_att_conv_head_0/bn/beta/Assign2layerfilter1PL_self_att_conv_head_0/bn/beta/read:02.layerfilter1PL_self_att_conv_head_0/bn/Const:08
Î
.layerfilter1PL_self_att_conv_head_0/bn/gamma:03layerfilter1PL_self_att_conv_head_0/bn/gamma/Assign3layerfilter1PL_self_att_conv_head_0/bn/gamma/read:020layerfilter1PL_self_att_conv_head_0/bn/Const_1:08

xlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assign}layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

zlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignlayerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1PL_self_att_conv_head_0/bn/layerfilter1PL_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
ã
-layerfilter1PL_neib_att_conv_head_0/weights:02layerfilter1PL_neib_att_conv_head_0/weights/Assign2layerfilter1PL_neib_att_conv_head_0/weights/read:02Hlayerfilter1PL_neib_att_conv_head_0/weights/Initializer/random_uniform:08
Ö
,layerfilter1PL_neib_att_conv_head_0/biases:01layerfilter1PL_neib_att_conv_head_0/biases/Assign1layerfilter1PL_neib_att_conv_head_0/biases/read:02>layerfilter1PL_neib_att_conv_head_0/biases/Initializer/Const:08
É
-layerfilter1PL_neib_att_conv_head_0/bn/beta:02layerfilter1PL_neib_att_conv_head_0/bn/beta/Assign2layerfilter1PL_neib_att_conv_head_0/bn/beta/read:02.layerfilter1PL_neib_att_conv_head_0/bn/Const:08
Î
.layerfilter1PL_neib_att_conv_head_0/bn/gamma:03layerfilter1PL_neib_att_conv_head_0/bn/gamma/Assign3layerfilter1PL_neib_att_conv_head_0/bn/gamma/read:020layerfilter1PL_neib_att_conv_head_0/bn/Const_1:08

xlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assign}layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

zlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignlayerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1PL_neib_att_conv_head_0/bn/layerfilter1PL_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
¢
layerfilter1PL/BiasAdd/biases:0$layerfilter1PL/BiasAdd/biases/Assign$layerfilter1PL/BiasAdd/biases/read:021layerfilter1PL/BiasAdd/biases/Initializer/zeros:08
w
gapnet10/weights:0gapnet10/weights/Assigngapnet10/weights/read:02-gapnet10/weights/Initializer/random_uniform:08
j
gapnet10/biases:0gapnet10/biases/Assigngapnet10/biases/read:02#gapnet10/biases/Initializer/Const:08
]
gapnet10/bn/beta:0gapnet10/bn/beta/Assigngapnet10/bn/beta/read:02gapnet10/bn/Const:08
b
gapnet10/bn/gamma:0gapnet10/bn/gamma/Assigngapnet10/bn/gamma/read:02gapnet10/bn/Const_1:08
¬
Bgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage:0Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/AssignGgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:02Tgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
Dgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage:0Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignIgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Vgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0

gapnet11PL/weights:0gapnet11PL/weights/Assigngapnet11PL/weights/read:02/gapnet11PL/weights/Initializer/random_uniform:08
r
gapnet11PL/biases:0gapnet11PL/biases/Assigngapnet11PL/biases/read:02%gapnet11PL/biases/Initializer/Const:08
e
gapnet11PL/bn/beta:0gapnet11PL/bn/beta/Assigngapnet11PL/bn/beta/read:02gapnet11PL/bn/Const:08
j
gapnet11PL/bn/gamma:0gapnet11PL/bn/gamma/Assigngapnet11PL/bn/gamma/read:02gapnet11PL/bn/Const_1:08
¼
Fgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage:0Kgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/AssignKgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/read:02Xgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
Ä
Hgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage:0Mgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignMgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Zgapnet11PL/bn/gapnet11PL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
k
aggPL/weights:0aggPL/weights/AssignaggPL/weights/read:02*aggPL/weights/Initializer/random_uniform:08
^
aggPL/biases:0aggPL/biases/AssignaggPL/biases/read:02 aggPL/biases/Initializer/Const:08
Q
aggPL/bn/beta:0aggPL/bn/beta/AssignaggPL/bn/beta/read:02aggPL/bn/Const:08
V
aggPL/bn/gamma:0aggPL/bn/gamma/AssignaggPL/bn/gamma/read:02aggPL/bn/Const_1:08

<aggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage:0AaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/AssignAaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/read:02NaggPL/bn/aggPL/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

>aggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage:0CaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignCaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02PaggPL/bn/aggPL/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
k
PLfc1/weights:0PLfc1/weights/AssignPLfc1/weights/read:02*PLfc1/weights/Initializer/random_uniform:08
^
PLfc1/biases:0PLfc1/biases/AssignPLfc1/biases/read:02 PLfc1/biases/Initializer/Const:08
Q
PLfc1/bn/beta:0PLfc1/bn/beta/AssignPLfc1/bn/beta/read:02PLfc1/bn/Const:08
V
PLfc1/bn/gamma:0PLfc1/bn/gamma/AssignPLfc1/bn/gamma/read:02PLfc1/bn/Const_1:08

<PLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage:0APLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/AssignAPLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/read:02NPLfc1/bn/PLfc1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

>PLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage:0CPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignCPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02PPLfc1/bn/PLfc1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
k
PLfc2/weights:0PLfc2/weights/AssignPLfc2/weights/read:02*PLfc2/weights/Initializer/random_uniform:08
^
PLfc2/biases:0PLfc2/biases/AssignPLfc2/biases/read:02 PLfc2/biases/Initializer/Const:08
Q
PLfc2/bn/beta:0PLfc2/bn/beta/AssignPLfc2/bn/beta/read:02PLfc2/bn/Const:08
V
PLfc2/bn/gamma:0PLfc2/bn/gamma/AssignPLfc2/bn/gamma/read:02PLfc2/bn/Const_1:08

<PLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage:0APLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/AssignAPLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/read:02NPLfc2/bn/PLfc2/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

>PLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage:0CPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignCPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02PPLfc2/bn/PLfc2/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
k
fc3PL/weights:0fc3PL/weights/Assignfc3PL/weights/read:02*fc3PL/weights/Initializer/random_uniform:08
^
fc3PL/biases:0fc3PL/biases/Assignfc3PL/biases/read:02 fc3PL/biases/Initializer/Const:08