óü-
µ))
,
Abs
x"T
y"T"
Ttype:

2	
:
Add
x"T
y"T
z"T"
Ttype:
2	

ArgMax

input"T
	dimension"Tidx
output"output_type" 
Ttype:
2	"
Tidxtype0:
2	"
output_typetype0	:
2	
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
¼
AvgPool

value"T
output"T"
ksize	list(int)(0"
strides	list(int)(0""
paddingstring:
SAMEVALID"-
data_formatstringNHWC:
NHWCNCHW"
Ttype:
2
h
BatchMatMul
x"T
y"T
output"T"
Ttype:
	2"
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
ë
Conv2D

input"T
filter"T
output"T"
Ttype:
2"
strides	list(int)"
use_cudnn_on_gpubool(""
paddingstring:
SAMEVALID"-
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
2	
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
,
Floor
x"T
y"T"
Ttype:
2

Gather
params"Tparams
indices"Tindices
output"Tparams"
validate_indicesbool("
Tparamstype"
Tindicestype:
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
;
Maximum
x"T
y"T
z"T"
Ttype:

2	
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
2	
.
Neg
x"T
y"T"
Ttype:

2	
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
2	
\
	RefSwitch
data"T
pred

output_false"T
output_true"T"	
Ttype
D
Relu
features"T
activations"T"
Ttype:
2	
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

2
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
a
Slice

input"T
begin"Index
size"Index
output"T"	
Ttype"
Indextype:
2	
9
Softmax
logits"T
softmax"T"
Ttype:
2

#SparseSoftmaxCrossEntropyWithLogits
features"T
labels"Tlabels	
loss"T
backprop"T"
Ttype:
2"
Tlabelstype0	:
2	
1
Square
x"T
y"T"
Ttype:

2	
G
SquaredDifference
x"T
y"T
z"T"
Ttype:

2	
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
2	
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
shared_namestring "serve*1.6.02unknownñ$
d
PlaceholderPlaceholder*
dtype0*
shape:2*"
_output_shapes
:2
^
Placeholder_1Placeholder*
dtype0*
shape
:2*
_output_shapes

:2
^
Placeholder_2Placeholder*
dtype0*
shape
:*
_output_shapes

:
X
Variable/initial_valueConst*
value	B : *
dtype0*
_output_shapes
: 
l
Variable
VariableV2*
shape: *
dtype0*
	container *
shared_name *
_output_shapes
: 
¢
Variable/AssignAssignVariableVariable/initial_value*
T0*
validate_shape(*
use_locking(*
_class
loc:@Variable*
_output_shapes
: 
a
Variable/readIdentityVariable*
T0*
_class
loc:@Variable*
_output_shapes
: 
N
Placeholder_3Placeholder*
dtype0
*
shape: *
_output_shapes
: 
^
Placeholder_4Placeholder*
dtype0*
shape
:2*
_output_shapes

:2
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
T0*

Tdim0*&
_output_shapes
:2
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
strided_sliceStridedSlicePlaceholderstrided_slice/stackstrided_slice/stack_1strided_slice/stack_2*
T0*
Index0*

begin_mask*
end_mask*
ellipsis_mask *
new_axis_mask *
shrink_axis_mask *"
_output_shapes
:2
^
SqueezeSqueezestrided_slice*
T0*
squeeze_dims
 *
_output_shapes

:2
R
ExpandDims_1/dimConst*
value	B : *
dtype0*
_output_shapes
: 
n
ExpandDims_1
ExpandDimsSqueezeExpandDims_1/dim*
T0*

Tdim0*"
_output_shapes
:2
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
strided_slice_1/stack_2Const*!
valueB"         *
dtype0*
_output_shapes
:

strided_slice_1StridedSliceExpandDims_1strided_slice_1/stackstrided_slice_1/stack_1strided_slice_1/stack_2*
T0*
Index0*

begin_mask*
end_mask*
ellipsis_mask *
new_axis_mask *
shrink_axis_mask*
_output_shapes

:2
[
ExpandDims_2/dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
v
ExpandDims_2
ExpandDimsstrided_slice_1ExpandDims_2/dim*
T0*

Tdim0*"
_output_shapes
:2
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
strided_slice_2/stack_2Const*
valueB"      *
dtype0*
_output_shapes
:

strided_slice_2StridedSliceExpandDims_1strided_slice_2/stackstrided_slice_2/stack_1strided_slice_2/stack_2*
T0*
Index0*

begin_mask*
end_mask*
ellipsis_mask *
new_axis_mask *
shrink_axis_mask *"
_output_shapes
:2
L
Equal/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
R
EqualEqualExpandDims_2Equal/y*
T0*"
_output_shapes
:2
d
ones_like/ShapeConst*!
valueB"   2      *
dtype0*
_output_shapes
:
T
ones_like/ConstConst*
valueB
 *  ?*
dtype0*
_output_shapes
: 
r
	ones_likeFillones_like/Shapeones_like/Const*
T0*

index_type0*"
_output_shapes
:2
Z
ShapeConst*!
valueB"   2      *
dtype0*
_output_shapes
:
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
:2
U
SelectSelectEqual	ones_likeFill*
T0*"
_output_shapes
:2
J
mul/xConst*
valueB
 *  zD*
dtype0*
_output_shapes
: 
F
mulMulmul/xSelect*
T0*"
_output_shapes
:2
c
transpose/permConst*!
valueB"          *
dtype0*
_output_shapes
:
e
	transpose	Transposemultranspose/perm*
T0*
Tperm0*"
_output_shapes
:2
j
strided_slice_3/stackConst*!
valueB"            *
dtype0*
_output_shapes
:
l
strided_slice_3/stack_1Const*!
valueB"           *
dtype0*
_output_shapes
:
l
strided_slice_3/stack_2Const*!
valueB"         *
dtype0*
_output_shapes
:

strided_slice_3StridedSliceExpandDims_1strided_slice_3/stackstrided_slice_3/stack_1strided_slice_3/stack_2*
T0*
Index0*

begin_mask*
end_mask*
ellipsis_mask *
new_axis_mask *
shrink_axis_mask *"
_output_shapes
:2
e
transpose_1/permConst*!
valueB"          *
dtype0*
_output_shapes
:
u
transpose_1	Transposestrided_slice_3transpose_1/perm*
T0*
Tperm0*"
_output_shapes
:2
j
strided_slice_4/stackConst*!
valueB"           *
dtype0*
_output_shapes
:
l
strided_slice_4/stack_1Const*!
valueB"            *
dtype0*
_output_shapes
:
l
strided_slice_4/stack_2Const*!
valueB"         *
dtype0*
_output_shapes
:

strided_slice_4StridedSlicetranspose_1strided_slice_4/stackstrided_slice_4/stack_1strided_slice_4/stack_2*
T0*
Index0*

begin_mask*
end_mask*
ellipsis_mask *
new_axis_mask *
shrink_axis_mask *"
_output_shapes
:2
c
Tile/multiplesConst*!
valueB"   2      *
dtype0*
_output_shapes
:
l
TileTilestrided_slice_4Tile/multiples*
T0*

Tmultiples0*"
_output_shapes
:22
e
transpose_2/permConst*!
valueB"          *
dtype0*
_output_shapes
:
j
transpose_2	TransposeTiletranspose_2/perm*
T0*
Tperm0*"
_output_shapes
:22
J
subSubTiletranspose_2*
T0*"
_output_shapes
:22
<
AbsAbssub*
T0*"
_output_shapes
:22
>
Abs_1AbsAbs*
T0*"
_output_shapes
:22
S
GreaterEqual/yConst*
valueB
 *ÛÉ@*
dtype0*
_output_shapes
: 
`
GreaterEqualGreaterEqualAbs_1GreaterEqual/y*
T0*"
_output_shapes
:22
L
mul_1/xConst*
valueB
 *ÛIA*
dtype0*
_output_shapes
: 
G
mul_1Mulmul_1/xAbs*
T0*"
_output_shapes
:22
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
:22
C
sub_2SubAbsAbs*
T0*"
_output_shapes
:22
[
Select_1SelectGreaterEqualsub_1sub_2*
T0*"
_output_shapes
:22
z
MatMulBatchMatMulstrided_slice_3transpose_1*
T0*
adj_x( *
adj_y( *"
_output_shapes
:22
L
mul_2/xConst*
valueB
 *   À*
dtype0*
_output_shapes
: 
J
mul_2Mulmul_2/xMatMul*
T0*"
_output_shapes
:22
N
SquareSquarestrided_slice_3*
T0*"
_output_shapes
:2
`
Sum/reduction_indicesConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
s
SumSumSquareSum/reduction_indices*
	keep_dims(*
T0*

Tidx0*"
_output_shapes
:2
e
transpose_3/permConst*!
valueB"          *
dtype0*
_output_shapes
:
i
transpose_3	TransposeSumtranspose_3/perm*
T0*
Tperm0*"
_output_shapes
:2
C
addAddSummul_2*
T0*"
_output_shapes
:22
K
add_1Addaddtranspose_3*
T0*"
_output_shapes
:22
J
add_2Addadd_1Select_1*
T0*"
_output_shapes
:22
E
add_3Addadd_2mul*
T0*"
_output_shapes
:22
K
add_4Addadd_3	transpose*
T0*"
_output_shapes
:22
>
NegNegadd_4*
T0*"
_output_shapes
:22
J
TopKV2/kConst*
value	B :
*
dtype0*
_output_shapes
: 
h
TopKV2TopKV2NegTopKV2/k*
sorted(*
T0*0
_output_shapes
:2
:2

^
	Squeeze_1SqueezePlaceholder*
T0*
squeeze_dims
 *
_output_shapes

:2
R
ExpandDims_3/dimConst*
value	B : *
dtype0*
_output_shapes
: 
p
ExpandDims_3
ExpandDims	Squeeze_1ExpandDims_3/dim*
T0*

Tdim0*"
_output_shapes
:2
[
ExpandDims_4/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
w
ExpandDims_4
ExpandDimsExpandDims_3ExpandDims_4/dim*
T0*

Tdim0*&
_output_shapes
:2
Ý
Hlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"             *
dtype0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*
_output_shapes
:
Ç
Flayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *øKÆ¾*
dtype0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*
_output_shapes
: 
Ç
Flayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *øKÆ>*
dtype0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*
_output_shapes
: 
¾
Playerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformHlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 
º
Flayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/subSubFlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/maxFlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/min*
T0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*
_output_shapes
: 
Ô
Flayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/mulMulPlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformFlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/sub*
T0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 
Æ
Blayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniformAddFlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/mulFlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform/min*
T0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 
ö
'layerfilter0_newfea_conv_head_0/weights
VariableV2"/device:CPU:0*
shape: *
dtype0*
	container *
shared_name *:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 
Ê
.layerfilter0_newfea_conv_head_0/weights/AssignAssign'layerfilter0_newfea_conv_head_0/weightsBlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 
Ý
,layerfilter0_newfea_conv_head_0/weights/readIdentity'layerfilter0_newfea_conv_head_0/weights"/device:CPU:0*
T0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 

&layerfilter0_newfea_conv_head_0/L2LossL2Loss,layerfilter0_newfea_conv_head_0/weights/read*
T0*
_output_shapes
: 
r
-layerfilter0_newfea_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
ª
+layerfilter0_newfea_conv_head_0/weight_lossMul&layerfilter0_newfea_conv_head_0/L2Loss-layerfilter0_newfea_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 

&layerfilter0_newfea_conv_head_0/Conv2DConv2DExpandDims_4,layerfilter0_newfea_conv_head_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2 
u
(layerfilter0_newfea_conv_head_0/bn/ConstConst*
valueB *    *
dtype0*
_output_shapes
: 

'layerfilter0_newfea_conv_head_0/bn/beta
VariableV2*
shape: *
dtype0*
	container *
shared_name *
_output_shapes
: 

.layerfilter0_newfea_conv_head_0/bn/beta/AssignAssign'layerfilter0_newfea_conv_head_0/bn/beta(layerfilter0_newfea_conv_head_0/bn/Const*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/bn/beta*
_output_shapes
: 
Â
,layerfilter0_newfea_conv_head_0/bn/beta/readIdentity'layerfilter0_newfea_conv_head_0/bn/beta*
T0*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/bn/beta*
_output_shapes
: 
w
*layerfilter0_newfea_conv_head_0/bn/Const_1Const*
valueB *  ?*
dtype0*
_output_shapes
: 

(layerfilter0_newfea_conv_head_0/bn/gamma
VariableV2*
shape: *
dtype0*
	container *
shared_name *
_output_shapes
: 

/layerfilter0_newfea_conv_head_0/bn/gamma/AssignAssign(layerfilter0_newfea_conv_head_0/bn/gamma*layerfilter0_newfea_conv_head_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_newfea_conv_head_0/bn/gamma*
_output_shapes
: 
Å
-layerfilter0_newfea_conv_head_0/bn/gamma/readIdentity(layerfilter0_newfea_conv_head_0/bn/gamma*
T0*;
_class1
/-loc:@layerfilter0_newfea_conv_head_0/bn/gamma*
_output_shapes
: 

Alayerfilter0_newfea_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ð
/layerfilter0_newfea_conv_head_0/bn/moments/meanMean&layerfilter0_newfea_conv_head_0/Conv2DAlayerfilter0_newfea_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
: 
©
7layerfilter0_newfea_conv_head_0/bn/moments/StopGradientStopGradient/layerfilter0_newfea_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
: 
ã
<layerfilter0_newfea_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference&layerfilter0_newfea_conv_head_0/Conv2D7layerfilter0_newfea_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2 

Elayerfilter0_newfea_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

3layerfilter0_newfea_conv_head_0/bn/moments/varianceMean<layerfilter0_newfea_conv_head_0/bn/moments/SquaredDifferenceElayerfilter0_newfea_conv_head_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
: 
¬
2layerfilter0_newfea_conv_head_0/bn/moments/SqueezeSqueeze/layerfilter0_newfea_conv_head_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
: 
²
4layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1Squeeze3layerfilter0_newfea_conv_head_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
: 
y
.layerfilter0_newfea_conv_head_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

0layerfilter0_newfea_conv_head_0/bn/cond/switch_tIdentity0layerfilter0_newfea_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

0layerfilter0_newfea_conv_head_0/bn/cond/switch_fIdentity.layerfilter0_newfea_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
k
/layerfilter0_newfea_conv_head_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ß
layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB: *
dtype0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ð
layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
à
layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Þ
nlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape: *
dtype0*
	container *
shared_name *
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ä
ulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

slayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentitynlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ã
layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB: *
dtype0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ô
layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
è
layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
â
playerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape: *
dtype0*
	container *
shared_name *
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ì
wlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

ulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¾
Flayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst1^layerfilter0_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ò
Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst1^layerfilter0_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
 
Tlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubVlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xFlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ê
Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Sub_layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1alayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¼
]layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchslayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read/layerfilter0_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
: : 
À
_layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch2layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/layerfilter0_newfea_conv_head_0/bn/cond/pred_id*
T0*E
_class;
97loc:@layerfilter0_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
: : 
²
Tlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulVlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Tlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ê
Playerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubYlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Tlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
´
Wlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/layerfilter0_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
: : 
Ö
Xlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst1^layerfilter0_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¦
Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubXlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xFlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ò
Xlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subalayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1clayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Â
_layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read/layerfilter0_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
: : 
Æ
alayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch4layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/layerfilter0_newfea_conv_head_0/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
: : 
º
Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulXlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ò
Rlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub[layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
º
Ylayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/layerfilter0_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
: : 
£
@layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverageNoOpQ^layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgS^layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_11^layerfilter0_newfea_conv_head_0/bn/cond/switch_t
¡
:layerfilter0_newfea_conv_head_0/bn/cond/control_dependencyIdentity0layerfilter0_newfea_conv_head_0/bn/cond/switch_tA^layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*C
_class9
75loc:@layerfilter0_newfea_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
g
,layerfilter0_newfea_conv_head_0/bn/cond/NoOpNoOp1^layerfilter0_newfea_conv_head_0/bn/cond/switch_f

<layerfilter0_newfea_conv_head_0/bn/cond/control_dependency_1Identity0layerfilter0_newfea_conv_head_0/bn/cond/switch_f-^layerfilter0_newfea_conv_head_0/bn/cond/NoOp*
T0
*C
_class9
75loc:@layerfilter0_newfea_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
Ü
-layerfilter0_newfea_conv_head_0/bn/cond/MergeMerge<layerfilter0_newfea_conv_head_0/bn/cond/control_dependency_1:layerfilter0_newfea_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
{
0layerfilter0_newfea_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter0_newfea_conv_head_0/bn/cond_1/switch_tIdentity2layerfilter0_newfea_conv_head_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

2layerfilter0_newfea_conv_head_0/bn/cond_1/switch_fIdentity0layerfilter0_newfea_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
m
1layerfilter0_newfea_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ð
2layerfilter0_newfea_conv_head_0/bn/cond_1/IdentityIdentity;layerfilter0_newfea_conv_head_0/bn/cond_1/Identity/Switch:1.^layerfilter0_newfea_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
: 

9layerfilter0_newfea_conv_head_0/bn/cond_1/Identity/SwitchSwitch2layerfilter0_newfea_conv_head_0/bn/moments/Squeeze1layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id*
T0*E
_class;
97loc:@layerfilter0_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
: : 
Ô
4layerfilter0_newfea_conv_head_0/bn/cond_1/Identity_1Identity=layerfilter0_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1.^layerfilter0_newfea_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
: 
¢
;layerfilter0_newfea_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch4layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_11layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
: : 

2layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_1Switchslayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id*
T0*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
: : 

2layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_2Switchulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id*
T0*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
: : 
Ð
/layerfilter0_newfea_conv_head_0/bn/cond_1/MergeMerge2layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_12layerfilter0_newfea_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

: : 
Ô
1layerfilter0_newfea_conv_head_0/bn/cond_1/Merge_1Merge2layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_24layerfilter0_newfea_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

: : 
w
2layerfilter0_newfea_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
Ã
0layerfilter0_newfea_conv_head_0/bn/batchnorm/addAdd1layerfilter0_newfea_conv_head_0/bn/cond_1/Merge_12layerfilter0_newfea_conv_head_0/bn/batchnorm/add/y*
T0*
_output_shapes
: 

2layerfilter0_newfea_conv_head_0/bn/batchnorm/RsqrtRsqrt0layerfilter0_newfea_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
: 
¿
0layerfilter0_newfea_conv_head_0/bn/batchnorm/mulMul2layerfilter0_newfea_conv_head_0/bn/batchnorm/Rsqrt-layerfilter0_newfea_conv_head_0/bn/gamma/read*
T0*
_output_shapes
: 
Ä
2layerfilter0_newfea_conv_head_0/bn/batchnorm/mul_1Mul&layerfilter0_newfea_conv_head_0/Conv2D0layerfilter0_newfea_conv_head_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2 
Á
2layerfilter0_newfea_conv_head_0/bn/batchnorm/mul_2Mul/layerfilter0_newfea_conv_head_0/bn/cond_1/Merge0layerfilter0_newfea_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
: 
¾
0layerfilter0_newfea_conv_head_0/bn/batchnorm/subSub,layerfilter0_newfea_conv_head_0/bn/beta/read2layerfilter0_newfea_conv_head_0/bn/batchnorm/mul_2*
T0*
_output_shapes
: 
Ð
2layerfilter0_newfea_conv_head_0/bn/batchnorm/add_1Add2layerfilter0_newfea_conv_head_0/bn/batchnorm/mul_10layerfilter0_newfea_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2 

$layerfilter0_newfea_conv_head_0/ReluRelu2layerfilter0_newfea_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2 
_
	Squeeze_2SqueezeExpandDims_4*
T0*
squeeze_dims
 *
_output_shapes

:2
R
ExpandDims_5/dimConst*
value	B : *
dtype0*
_output_shapes
: 
p
ExpandDims_5
ExpandDims	Squeeze_2ExpandDims_5/dim*
T0*

Tdim0*"
_output_shapes
:2
M
range/startConst*
value	B : *
dtype0*
_output_shapes
: 
M
range/limitConst*
value	B :*
dtype0*
_output_shapes
: 
M
range/deltaConst*
value	B :*
dtype0*
_output_shapes
: 
]
rangeRangerange/startrange/limitrange/delta*

Tidx0*
_output_shapes
:
I
mul_3/yConst*
value	B :2*
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
valueB"ÿÿÿÿ   *
dtype0*
_output_shapes
:
j
	Reshape_1ReshapeExpandDims_5Reshape_1/shape*
T0*
Tshape0*
_output_shapes

:2
L
add_5AddTopKV2:1Reshape*
T0*"
_output_shapes
:2


GatherGather	Reshape_1add_5*
validate_indices(*
Tparams0*
Tindices0*&
_output_shapes
:2

i
Tile_1/multiplesConst*%
valueB"      
      *
dtype0*
_output_shapes
:
q
Tile_1TileExpandDims_4Tile_1/multiples*
T0*

Tmultiples0*&
_output_shapes
:2

M
sub_3SubTile_1Gather*
T0*&
_output_shapes
:2

Ë
?layerfilter0_edgefea_0/weights/Initializer/random_uniform/shapeConst*%
valueB"             *
dtype0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*
_output_shapes
:
µ
=layerfilter0_edgefea_0/weights/Initializer/random_uniform/minConst*
valueB
 *øKÆ¾*
dtype0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*
_output_shapes
: 
µ
=layerfilter0_edgefea_0/weights/Initializer/random_uniform/maxConst*
valueB
 *øKÆ>*
dtype0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*
_output_shapes
: 
£
Glayerfilter0_edgefea_0/weights/Initializer/random_uniform/RandomUniformRandomUniform?layerfilter0_edgefea_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 

=layerfilter0_edgefea_0/weights/Initializer/random_uniform/subSub=layerfilter0_edgefea_0/weights/Initializer/random_uniform/max=layerfilter0_edgefea_0/weights/Initializer/random_uniform/min*
T0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*
_output_shapes
: 
°
=layerfilter0_edgefea_0/weights/Initializer/random_uniform/mulMulGlayerfilter0_edgefea_0/weights/Initializer/random_uniform/RandomUniform=layerfilter0_edgefea_0/weights/Initializer/random_uniform/sub*
T0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 
¢
9layerfilter0_edgefea_0/weights/Initializer/random_uniformAdd=layerfilter0_edgefea_0/weights/Initializer/random_uniform/mul=layerfilter0_edgefea_0/weights/Initializer/random_uniform/min*
T0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 
ä
layerfilter0_edgefea_0/weights
VariableV2"/device:CPU:0*
shape: *
dtype0*
	container *
shared_name *1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 
¦
%layerfilter0_edgefea_0/weights/AssignAssignlayerfilter0_edgefea_0/weights9layerfilter0_edgefea_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 
Â
#layerfilter0_edgefea_0/weights/readIdentitylayerfilter0_edgefea_0/weights"/device:CPU:0*
T0*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 
m
layerfilter0_edgefea_0/L2LossL2Loss#layerfilter0_edgefea_0/weights/read*
T0*
_output_shapes
: 
i
$layerfilter0_edgefea_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 

"layerfilter0_edgefea_0/weight_lossMullayerfilter0_edgefea_0/L2Loss$layerfilter0_edgefea_0/weight_loss/y*
T0*
_output_shapes
: 
ó
layerfilter0_edgefea_0/Conv2DConv2Dsub_3#layerfilter0_edgefea_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2
 
®
/layerfilter0_edgefea_0/biases/Initializer/ConstConst*
valueB *    *
dtype0*0
_class&
$"loc:@layerfilter0_edgefea_0/biases*
_output_shapes
: 
Ê
layerfilter0_edgefea_0/biases
VariableV2"/device:CPU:0*
shape: *
dtype0*
	container *
shared_name *0
_class&
$"loc:@layerfilter0_edgefea_0/biases*
_output_shapes
: 

$layerfilter0_edgefea_0/biases/AssignAssignlayerfilter0_edgefea_0/biases/layerfilter0_edgefea_0/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter0_edgefea_0/biases*
_output_shapes
: 
³
"layerfilter0_edgefea_0/biases/readIdentitylayerfilter0_edgefea_0/biases"/device:CPU:0*
T0*0
_class&
$"loc:@layerfilter0_edgefea_0/biases*
_output_shapes
: 
´
layerfilter0_edgefea_0/BiasAddBiasAddlayerfilter0_edgefea_0/Conv2D"layerfilter0_edgefea_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
 
l
layerfilter0_edgefea_0/bn/ConstConst*
valueB *    *
dtype0*
_output_shapes
: 

layerfilter0_edgefea_0/bn/beta
VariableV2*
shape: *
dtype0*
	container *
shared_name *
_output_shapes
: 
ñ
%layerfilter0_edgefea_0/bn/beta/AssignAssignlayerfilter0_edgefea_0/bn/betalayerfilter0_edgefea_0/bn/Const*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter0_edgefea_0/bn/beta*
_output_shapes
: 
§
#layerfilter0_edgefea_0/bn/beta/readIdentitylayerfilter0_edgefea_0/bn/beta*
T0*1
_class'
%#loc:@layerfilter0_edgefea_0/bn/beta*
_output_shapes
: 
n
!layerfilter0_edgefea_0/bn/Const_1Const*
valueB *  ?*
dtype0*
_output_shapes
: 

layerfilter0_edgefea_0/bn/gamma
VariableV2*
shape: *
dtype0*
	container *
shared_name *
_output_shapes
: 
ö
&layerfilter0_edgefea_0/bn/gamma/AssignAssignlayerfilter0_edgefea_0/bn/gamma!layerfilter0_edgefea_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter0_edgefea_0/bn/gamma*
_output_shapes
: 
ª
$layerfilter0_edgefea_0/bn/gamma/readIdentitylayerfilter0_edgefea_0/bn/gamma*
T0*2
_class(
&$loc:@layerfilter0_edgefea_0/bn/gamma*
_output_shapes
: 

8layerfilter0_edgefea_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ö
&layerfilter0_edgefea_0/bn/moments/meanMeanlayerfilter0_edgefea_0/BiasAdd8layerfilter0_edgefea_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
: 

.layerfilter0_edgefea_0/bn/moments/StopGradientStopGradient&layerfilter0_edgefea_0/bn/moments/mean*
T0*&
_output_shapes
: 
É
3layerfilter0_edgefea_0/bn/moments/SquaredDifferenceSquaredDifferencelayerfilter0_edgefea_0/BiasAdd.layerfilter0_edgefea_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2
 

<layerfilter0_edgefea_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ó
*layerfilter0_edgefea_0/bn/moments/varianceMean3layerfilter0_edgefea_0/bn/moments/SquaredDifference<layerfilter0_edgefea_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
: 

)layerfilter0_edgefea_0/bn/moments/SqueezeSqueeze&layerfilter0_edgefea_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
: 
 
+layerfilter0_edgefea_0/bn/moments/Squeeze_1Squeeze*layerfilter0_edgefea_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
: 
p
%layerfilter0_edgefea_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
}
'layerfilter0_edgefea_0/bn/cond/switch_tIdentity'layerfilter0_edgefea_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 
{
'layerfilter0_edgefea_0/bn/cond/switch_fIdentity%layerfilter0_edgefea_0/bn/cond/Switch*
T0
*
_output_shapes
: 
b
&layerfilter0_edgefea_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
¹
~layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB: *
dtype0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ª
tlayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

nlayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFill~layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensortlayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¹
\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape: *
dtype0*
	container *
shared_name *o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ú
clayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragenlayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
á
alayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¾
layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB: *
dtype0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
®
vlayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

playerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorvlayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
½
^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape: *
dtype0*
	container *
shared_name *q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

elayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssign^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageplayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ç
clayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentity^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¬
=layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/decayConst(^layerfilter0_edgefea_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
­
Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst(^layerfilter0_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ò
Klayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubMlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x=layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
T0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubVlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Xlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

Tlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchalayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read&layerfilter0_edgefea_0/bn/cond/pred_id*
T0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
: : 

Vlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch)layerfilter0_edgefea_0/bn/moments/Squeeze&layerfilter0_edgefea_0/bn/cond/pred_id*
T0*<
_class2
0.loc:@layerfilter0_edgefea_0/bn/moments/Squeeze* 
_output_shapes
: : 

Klayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulMlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Klayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

Glayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubPlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Klayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ý
Nlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage&layerfilter0_edgefea_0/bn/cond/pred_id*
T0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
: : 
±
Olayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst(^layerfilter0_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubOlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x=layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
T0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¤
Olayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubXlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Zlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

Vlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchclayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read&layerfilter0_edgefea_0/bn/cond/pred_id*
T0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
: : 
¢
Xlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch+layerfilter0_edgefea_0/bn/moments/Squeeze_1&layerfilter0_edgefea_0/bn/cond/pred_id*
T0*>
_class4
20loc:@layerfilter0_edgefea_0/bn/moments/Squeeze_1* 
_output_shapes
: : 

Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulOlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¤
Ilayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubRlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

Playerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitch^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage&layerfilter0_edgefea_0/bn/cond/pred_id*
T0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
: : 
ÿ
7layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverageNoOpH^layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgJ^layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1(^layerfilter0_edgefea_0/bn/cond/switch_t
ý
1layerfilter0_edgefea_0/bn/cond/control_dependencyIdentity'layerfilter0_edgefea_0/bn/cond/switch_t8^layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage*
T0
*:
_class0
.,loc:@layerfilter0_edgefea_0/bn/cond/switch_t*
_output_shapes
: 
U
#layerfilter0_edgefea_0/bn/cond/NoOpNoOp(^layerfilter0_edgefea_0/bn/cond/switch_f
ë
3layerfilter0_edgefea_0/bn/cond/control_dependency_1Identity'layerfilter0_edgefea_0/bn/cond/switch_f$^layerfilter0_edgefea_0/bn/cond/NoOp*
T0
*:
_class0
.,loc:@layerfilter0_edgefea_0/bn/cond/switch_f*
_output_shapes
: 
Á
$layerfilter0_edgefea_0/bn/cond/MergeMerge3layerfilter0_edgefea_0/bn/cond/control_dependency_11layerfilter0_edgefea_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
r
'layerfilter0_edgefea_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

)layerfilter0_edgefea_0/bn/cond_1/switch_tIdentity)layerfilter0_edgefea_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

)layerfilter0_edgefea_0/bn/cond_1/switch_fIdentity'layerfilter0_edgefea_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
d
(layerfilter0_edgefea_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
µ
)layerfilter0_edgefea_0/bn/cond_1/IdentityIdentity2layerfilter0_edgefea_0/bn/cond_1/Identity/Switch:1%^layerfilter0_edgefea_0/bn/cond/Merge*
T0*
_output_shapes
: 
ø
0layerfilter0_edgefea_0/bn/cond_1/Identity/SwitchSwitch)layerfilter0_edgefea_0/bn/moments/Squeeze(layerfilter0_edgefea_0/bn/cond_1/pred_id*
T0*<
_class2
0.loc:@layerfilter0_edgefea_0/bn/moments/Squeeze* 
_output_shapes
: : 
¹
+layerfilter0_edgefea_0/bn/cond_1/Identity_1Identity4layerfilter0_edgefea_0/bn/cond_1/Identity_1/Switch:1%^layerfilter0_edgefea_0/bn/cond/Merge*
T0*
_output_shapes
: 
þ
2layerfilter0_edgefea_0/bn/cond_1/Identity_1/SwitchSwitch+layerfilter0_edgefea_0/bn/moments/Squeeze_1(layerfilter0_edgefea_0/bn/cond_1/pred_id*
T0*>
_class4
20loc:@layerfilter0_edgefea_0/bn/moments/Squeeze_1* 
_output_shapes
: : 
Ü
)layerfilter0_edgefea_0/bn/cond_1/Switch_1Switchalayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read(layerfilter0_edgefea_0/bn/cond_1/pred_id*
T0*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
: : 
à
)layerfilter0_edgefea_0/bn/cond_1/Switch_2Switchclayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read(layerfilter0_edgefea_0/bn/cond_1/pred_id*
T0*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
: : 
µ
&layerfilter0_edgefea_0/bn/cond_1/MergeMerge)layerfilter0_edgefea_0/bn/cond_1/Switch_1)layerfilter0_edgefea_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

: : 
¹
(layerfilter0_edgefea_0/bn/cond_1/Merge_1Merge)layerfilter0_edgefea_0/bn/cond_1/Switch_2+layerfilter0_edgefea_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

: : 
n
)layerfilter0_edgefea_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
¨
'layerfilter0_edgefea_0/bn/batchnorm/addAdd(layerfilter0_edgefea_0/bn/cond_1/Merge_1)layerfilter0_edgefea_0/bn/batchnorm/add/y*
T0*
_output_shapes
: 

)layerfilter0_edgefea_0/bn/batchnorm/RsqrtRsqrt'layerfilter0_edgefea_0/bn/batchnorm/add*
T0*
_output_shapes
: 
¤
'layerfilter0_edgefea_0/bn/batchnorm/mulMul)layerfilter0_edgefea_0/bn/batchnorm/Rsqrt$layerfilter0_edgefea_0/bn/gamma/read*
T0*
_output_shapes
: 
ª
)layerfilter0_edgefea_0/bn/batchnorm/mul_1Mullayerfilter0_edgefea_0/BiasAdd'layerfilter0_edgefea_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2
 
¦
)layerfilter0_edgefea_0/bn/batchnorm/mul_2Mul&layerfilter0_edgefea_0/bn/cond_1/Merge'layerfilter0_edgefea_0/bn/batchnorm/mul*
T0*
_output_shapes
: 
£
'layerfilter0_edgefea_0/bn/batchnorm/subSub#layerfilter0_edgefea_0/bn/beta/read)layerfilter0_edgefea_0/bn/batchnorm/mul_2*
T0*
_output_shapes
: 
µ
)layerfilter0_edgefea_0/bn/batchnorm/add_1Add)layerfilter0_edgefea_0/bn/batchnorm/mul_1'layerfilter0_edgefea_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2
 

layerfilter0_edgefea_0/ReluRelu)layerfilter0_edgefea_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2
 
á
Jlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"             *
dtype0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*
_output_shapes
:
Ë
Hlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *JQÚ¾*
dtype0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*
_output_shapes
: 
Ë
Hlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *JQÚ>*
dtype0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*
_output_shapes
: 
Ä
Rlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 
Â
Hlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/subSubHlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/maxHlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*
_output_shapes
: 
Ü
Hlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/mulMulRlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformHlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/sub*
T0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 
Î
Dlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniformAddHlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/mulHlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 
ú
)layerfilter0_self_att_conv_head_0/weights
VariableV2"/device:CPU:0*
shape: *
dtype0*
	container *
shared_name *<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 
Ò
0layerfilter0_self_att_conv_head_0/weights/AssignAssign)layerfilter0_self_att_conv_head_0/weightsDlayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 
ã
.layerfilter0_self_att_conv_head_0/weights/readIdentity)layerfilter0_self_att_conv_head_0/weights"/device:CPU:0*
T0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 

(layerfilter0_self_att_conv_head_0/L2LossL2Loss.layerfilter0_self_att_conv_head_0/weights/read*
T0*
_output_shapes
: 
t
/layerfilter0_self_att_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter0_self_att_conv_head_0/weight_lossMul(layerfilter0_self_att_conv_head_0/L2Loss/layerfilter0_self_att_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 
¨
(layerfilter0_self_att_conv_head_0/Conv2DConv2D$layerfilter0_newfea_conv_head_0/Relu.layerfilter0_self_att_conv_head_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2
Ä
:layerfilter0_self_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *
dtype0*;
_class1
/-loc:@layerfilter0_self_att_conv_head_0/biases*
_output_shapes
:
à
(layerfilter0_self_att_conv_head_0/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *;
_class1
/-loc:@layerfilter0_self_att_conv_head_0/biases*
_output_shapes
:
¹
/layerfilter0_self_att_conv_head_0/biases/AssignAssign(layerfilter0_self_att_conv_head_0/biases:layerfilter0_self_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_self_att_conv_head_0/biases*
_output_shapes
:
Ô
-layerfilter0_self_att_conv_head_0/biases/readIdentity(layerfilter0_self_att_conv_head_0/biases"/device:CPU:0*
T0*;
_class1
/-loc:@layerfilter0_self_att_conv_head_0/biases*
_output_shapes
:
Õ
)layerfilter0_self_att_conv_head_0/BiasAddBiasAdd(layerfilter0_self_att_conv_head_0/Conv2D-layerfilter0_self_att_conv_head_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
w
*layerfilter0_self_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

)layerfilter0_self_att_conv_head_0/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:

0layerfilter0_self_att_conv_head_0/bn/beta/AssignAssign)layerfilter0_self_att_conv_head_0/bn/beta*layerfilter0_self_att_conv_head_0/bn/Const*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/bn/beta*
_output_shapes
:
È
.layerfilter0_self_att_conv_head_0/bn/beta/readIdentity)layerfilter0_self_att_conv_head_0/bn/beta*
T0*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/bn/beta*
_output_shapes
:
y
,layerfilter0_self_att_conv_head_0/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

*layerfilter0_self_att_conv_head_0/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
¢
1layerfilter0_self_att_conv_head_0/bn/gamma/AssignAssign*layerfilter0_self_att_conv_head_0/bn/gamma,layerfilter0_self_att_conv_head_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter0_self_att_conv_head_0/bn/gamma*
_output_shapes
:
Ë
/layerfilter0_self_att_conv_head_0/bn/gamma/readIdentity*layerfilter0_self_att_conv_head_0/bn/gamma*
T0*=
_class3
1/loc:@layerfilter0_self_att_conv_head_0/bn/gamma*
_output_shapes
:

Clayerfilter0_self_att_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
÷
1layerfilter0_self_att_conv_head_0/bn/moments/meanMean)layerfilter0_self_att_conv_head_0/BiasAddClayerfilter0_self_att_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
­
9layerfilter0_self_att_conv_head_0/bn/moments/StopGradientStopGradient1layerfilter0_self_att_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
:
ê
>layerfilter0_self_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference)layerfilter0_self_att_conv_head_0/BiasAdd9layerfilter0_self_att_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2

Glayerfilter0_self_att_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

5layerfilter0_self_att_conv_head_0/bn/moments/varianceMean>layerfilter0_self_att_conv_head_0/bn/moments/SquaredDifferenceGlayerfilter0_self_att_conv_head_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
°
4layerfilter0_self_att_conv_head_0/bn/moments/SqueezeSqueeze1layerfilter0_self_att_conv_head_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:
¶
6layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1Squeeze5layerfilter0_self_att_conv_head_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:
{
0layerfilter0_self_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter0_self_att_conv_head_0/bn/cond/switch_tIdentity2layerfilter0_self_att_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter0_self_att_conv_head_0/bn/cond/switch_fIdentity0layerfilter0_self_att_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter0_self_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ç
layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ø
layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ð
layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
æ
rlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
ylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¤
wlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ë
layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ê
tlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
{layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ª
ylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Â
Hlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter0_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ú
Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter0_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ª
Vlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ô
Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
È
_layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter0_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
È
alayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter0_self_att_conv_head_0/bn/moments/Squeeze1layerfilter0_self_att_conv_head_0/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter0_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
¼
Vlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
Rlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
À
Ylayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter0_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Þ
Zlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter0_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
°
Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ü
Zlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Î
alayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter0_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Î
clayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_11layerfilter0_self_att_conv_head_0/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
Ä
Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
Tlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
[layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter0_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
«
Blayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpS^layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_13^layerfilter0_self_att_conv_head_0/bn/cond/switch_t
©
<layerfilter0_self_att_conv_head_0/bn/cond/control_dependencyIdentity2layerfilter0_self_att_conv_head_0/bn/cond/switch_tC^layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter0_self_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter0_self_att_conv_head_0/bn/cond/NoOpNoOp3^layerfilter0_self_att_conv_head_0/bn/cond/switch_f

>layerfilter0_self_att_conv_head_0/bn/cond/control_dependency_1Identity2layerfilter0_self_att_conv_head_0/bn/cond/switch_f/^layerfilter0_self_att_conv_head_0/bn/cond/NoOp*
T0
*E
_class;
97loc:@layerfilter0_self_att_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
â
/layerfilter0_self_att_conv_head_0/bn/cond/MergeMerge>layerfilter0_self_att_conv_head_0/bn/cond/control_dependency_1<layerfilter0_self_att_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter0_self_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

4layerfilter0_self_att_conv_head_0/bn/cond_1/switch_tIdentity4layerfilter0_self_att_conv_head_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

4layerfilter0_self_att_conv_head_0/bn/cond_1/switch_fIdentity2layerfilter0_self_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
o
3layerfilter0_self_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ö
4layerfilter0_self_att_conv_head_0/bn/cond_1/IdentityIdentity=layerfilter0_self_att_conv_head_0/bn/cond_1/Identity/Switch:10^layerfilter0_self_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
¤
;layerfilter0_self_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch4layerfilter0_self_att_conv_head_0/bn/moments/Squeeze3layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter0_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
Ú
6layerfilter0_self_att_conv_head_0/bn/cond_1/Identity_1Identity?layerfilter0_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:10^layerfilter0_self_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
ª
=layerfilter0_self_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch6layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_13layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::

4layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_1Switchwlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
£
4layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_2Switchylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
1layerfilter0_self_att_conv_head_0/bn/cond_1/MergeMerge4layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_14layerfilter0_self_att_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
Ú
3layerfilter0_self_att_conv_head_0/bn/cond_1/Merge_1Merge4layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_26layerfilter0_self_att_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
y
4layerfilter0_self_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
É
2layerfilter0_self_att_conv_head_0/bn/batchnorm/addAdd3layerfilter0_self_att_conv_head_0/bn/cond_1/Merge_14layerfilter0_self_att_conv_head_0/bn/batchnorm/add/y*
T0*
_output_shapes
:

4layerfilter0_self_att_conv_head_0/bn/batchnorm/RsqrtRsqrt2layerfilter0_self_att_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
:
Å
2layerfilter0_self_att_conv_head_0/bn/batchnorm/mulMul4layerfilter0_self_att_conv_head_0/bn/batchnorm/Rsqrt/layerfilter0_self_att_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:
Ë
4layerfilter0_self_att_conv_head_0/bn/batchnorm/mul_1Mul)layerfilter0_self_att_conv_head_0/BiasAdd2layerfilter0_self_att_conv_head_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2
Ç
4layerfilter0_self_att_conv_head_0/bn/batchnorm/mul_2Mul1layerfilter0_self_att_conv_head_0/bn/cond_1/Merge2layerfilter0_self_att_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:
Ä
2layerfilter0_self_att_conv_head_0/bn/batchnorm/subSub.layerfilter0_self_att_conv_head_0/bn/beta/read4layerfilter0_self_att_conv_head_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:
Ö
4layerfilter0_self_att_conv_head_0/bn/batchnorm/add_1Add4layerfilter0_self_att_conv_head_0/bn/batchnorm/mul_12layerfilter0_self_att_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2

&layerfilter0_self_att_conv_head_0/ReluRelu4layerfilter0_self_att_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2
á
Jlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"             *
dtype0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*
_output_shapes
:
Ë
Hlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *JQÚ¾*
dtype0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*
_output_shapes
: 
Ë
Hlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *JQÚ>*
dtype0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*
_output_shapes
: 
Ä
Rlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 
Â
Hlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/subSubHlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/maxHlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*
_output_shapes
: 
Ü
Hlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/mulMulRlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformHlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/sub*
T0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 
Î
Dlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniformAddHlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/mulHlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 
ú
)layerfilter0_neib_att_conv_head_0/weights
VariableV2"/device:CPU:0*
shape: *
dtype0*
	container *
shared_name *<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 
Ò
0layerfilter0_neib_att_conv_head_0/weights/AssignAssign)layerfilter0_neib_att_conv_head_0/weightsDlayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 
ã
.layerfilter0_neib_att_conv_head_0/weights/readIdentity)layerfilter0_neib_att_conv_head_0/weights"/device:CPU:0*
T0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 

(layerfilter0_neib_att_conv_head_0/L2LossL2Loss.layerfilter0_neib_att_conv_head_0/weights/read*
T0*
_output_shapes
: 
t
/layerfilter0_neib_att_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter0_neib_att_conv_head_0/weight_lossMul(layerfilter0_neib_att_conv_head_0/L2Loss/layerfilter0_neib_att_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 

(layerfilter0_neib_att_conv_head_0/Conv2DConv2Dlayerfilter0_edgefea_0/Relu.layerfilter0_neib_att_conv_head_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2

Ä
:layerfilter0_neib_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *
dtype0*;
_class1
/-loc:@layerfilter0_neib_att_conv_head_0/biases*
_output_shapes
:
à
(layerfilter0_neib_att_conv_head_0/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *;
_class1
/-loc:@layerfilter0_neib_att_conv_head_0/biases*
_output_shapes
:
¹
/layerfilter0_neib_att_conv_head_0/biases/AssignAssign(layerfilter0_neib_att_conv_head_0/biases:layerfilter0_neib_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_neib_att_conv_head_0/biases*
_output_shapes
:
Ô
-layerfilter0_neib_att_conv_head_0/biases/readIdentity(layerfilter0_neib_att_conv_head_0/biases"/device:CPU:0*
T0*;
_class1
/-loc:@layerfilter0_neib_att_conv_head_0/biases*
_output_shapes
:
Õ
)layerfilter0_neib_att_conv_head_0/BiasAddBiasAdd(layerfilter0_neib_att_conv_head_0/Conv2D-layerfilter0_neib_att_conv_head_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2

w
*layerfilter0_neib_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

)layerfilter0_neib_att_conv_head_0/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:

0layerfilter0_neib_att_conv_head_0/bn/beta/AssignAssign)layerfilter0_neib_att_conv_head_0/bn/beta*layerfilter0_neib_att_conv_head_0/bn/Const*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/bn/beta*
_output_shapes
:
È
.layerfilter0_neib_att_conv_head_0/bn/beta/readIdentity)layerfilter0_neib_att_conv_head_0/bn/beta*
T0*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/bn/beta*
_output_shapes
:
y
,layerfilter0_neib_att_conv_head_0/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

*layerfilter0_neib_att_conv_head_0/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
¢
1layerfilter0_neib_att_conv_head_0/bn/gamma/AssignAssign*layerfilter0_neib_att_conv_head_0/bn/gamma,layerfilter0_neib_att_conv_head_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter0_neib_att_conv_head_0/bn/gamma*
_output_shapes
:
Ë
/layerfilter0_neib_att_conv_head_0/bn/gamma/readIdentity*layerfilter0_neib_att_conv_head_0/bn/gamma*
T0*=
_class3
1/loc:@layerfilter0_neib_att_conv_head_0/bn/gamma*
_output_shapes
:

Clayerfilter0_neib_att_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
÷
1layerfilter0_neib_att_conv_head_0/bn/moments/meanMean)layerfilter0_neib_att_conv_head_0/BiasAddClayerfilter0_neib_att_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
­
9layerfilter0_neib_att_conv_head_0/bn/moments/StopGradientStopGradient1layerfilter0_neib_att_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
:
ê
>layerfilter0_neib_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference)layerfilter0_neib_att_conv_head_0/BiasAdd9layerfilter0_neib_att_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2


Glayerfilter0_neib_att_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

5layerfilter0_neib_att_conv_head_0/bn/moments/varianceMean>layerfilter0_neib_att_conv_head_0/bn/moments/SquaredDifferenceGlayerfilter0_neib_att_conv_head_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
°
4layerfilter0_neib_att_conv_head_0/bn/moments/SqueezeSqueeze1layerfilter0_neib_att_conv_head_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:
¶
6layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1Squeeze5layerfilter0_neib_att_conv_head_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:
{
0layerfilter0_neib_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter0_neib_att_conv_head_0/bn/cond/switch_tIdentity2layerfilter0_neib_att_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter0_neib_att_conv_head_0/bn/cond/switch_fIdentity0layerfilter0_neib_att_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter0_neib_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ç
layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ø
layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ð
layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
æ
rlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
ylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¤
wlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ë
layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ê
tlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
{layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ª
ylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Â
Hlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter0_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ú
Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter0_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ª
Vlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ô
Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
È
_layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter0_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
È
alayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze1layerfilter0_neib_att_conv_head_0/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
¼
Vlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
Rlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
À
Ylayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter0_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Þ
Zlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter0_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
°
Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ü
Zlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Î
alayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter0_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Î
clayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_11layerfilter0_neib_att_conv_head_0/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
Ä
Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
Tlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
[layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter0_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
«
Blayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpS^layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_13^layerfilter0_neib_att_conv_head_0/bn/cond/switch_t
©
<layerfilter0_neib_att_conv_head_0/bn/cond/control_dependencyIdentity2layerfilter0_neib_att_conv_head_0/bn/cond/switch_tC^layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter0_neib_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter0_neib_att_conv_head_0/bn/cond/NoOpNoOp3^layerfilter0_neib_att_conv_head_0/bn/cond/switch_f

>layerfilter0_neib_att_conv_head_0/bn/cond/control_dependency_1Identity2layerfilter0_neib_att_conv_head_0/bn/cond/switch_f/^layerfilter0_neib_att_conv_head_0/bn/cond/NoOp*
T0
*E
_class;
97loc:@layerfilter0_neib_att_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
â
/layerfilter0_neib_att_conv_head_0/bn/cond/MergeMerge>layerfilter0_neib_att_conv_head_0/bn/cond/control_dependency_1<layerfilter0_neib_att_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter0_neib_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

4layerfilter0_neib_att_conv_head_0/bn/cond_1/switch_tIdentity4layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

4layerfilter0_neib_att_conv_head_0/bn/cond_1/switch_fIdentity2layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
o
3layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ö
4layerfilter0_neib_att_conv_head_0/bn/cond_1/IdentityIdentity=layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity/Switch:10^layerfilter0_neib_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
¤
;layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch4layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze3layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
Ú
6layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity_1Identity?layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:10^layerfilter0_neib_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
ª
=layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch6layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_13layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::

4layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_1Switchwlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
£
4layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_2Switchylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
1layerfilter0_neib_att_conv_head_0/bn/cond_1/MergeMerge4layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_14layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
Ú
3layerfilter0_neib_att_conv_head_0/bn/cond_1/Merge_1Merge4layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_26layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
y
4layerfilter0_neib_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
É
2layerfilter0_neib_att_conv_head_0/bn/batchnorm/addAdd3layerfilter0_neib_att_conv_head_0/bn/cond_1/Merge_14layerfilter0_neib_att_conv_head_0/bn/batchnorm/add/y*
T0*
_output_shapes
:

4layerfilter0_neib_att_conv_head_0/bn/batchnorm/RsqrtRsqrt2layerfilter0_neib_att_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
:
Å
2layerfilter0_neib_att_conv_head_0/bn/batchnorm/mulMul4layerfilter0_neib_att_conv_head_0/bn/batchnorm/Rsqrt/layerfilter0_neib_att_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:
Ë
4layerfilter0_neib_att_conv_head_0/bn/batchnorm/mul_1Mul)layerfilter0_neib_att_conv_head_0/BiasAdd2layerfilter0_neib_att_conv_head_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2

Ç
4layerfilter0_neib_att_conv_head_0/bn/batchnorm/mul_2Mul1layerfilter0_neib_att_conv_head_0/bn/cond_1/Merge2layerfilter0_neib_att_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:
Ä
2layerfilter0_neib_att_conv_head_0/bn/batchnorm/subSub.layerfilter0_neib_att_conv_head_0/bn/beta/read4layerfilter0_neib_att_conv_head_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:
Ö
4layerfilter0_neib_att_conv_head_0/bn/batchnorm/add_1Add4layerfilter0_neib_att_conv_head_0/bn/batchnorm/mul_12layerfilter0_neib_att_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2


&layerfilter0_neib_att_conv_head_0/ReluRelu4layerfilter0_neib_att_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2


add_6Add&layerfilter0_self_att_conv_head_0/Relu&layerfilter0_neib_att_conv_head_0/Relu*
T0*&
_output_shapes
:2

i
transpose_4/permConst*%
valueB"             *
dtype0*
_output_shapes
:
o
transpose_4	Transposeadd_6transpose_4/perm*
T0*
Tperm0*&
_output_shapes
:2

T
LeakyRelu/alphaConst*
valueB
 *ÍÌL>*
dtype0*
_output_shapes
: 
c
LeakyRelu/mulMulLeakyRelu/alphatranspose_4*
T0*&
_output_shapes
:2

i
LeakyRelu/MaximumMaximumLeakyRelu/multranspose_4*
T0*&
_output_shapes
:2

`
Shape_1Const*%
valueB"   2      
   *
dtype0*
_output_shapes
:
F
RankConst*
value	B :*
dtype0*
_output_shapes
: 
`
Shape_2Const*%
valueB"   2      
   *
dtype0*
_output_shapes
:
G
Sub/yConst*
value	B :*
dtype0*
_output_shapes
: 
8
SubSubRankSub/y*
T0*
_output_shapes
: 
R
Slice/beginPackSub*
N*
T0*

axis *
_output_shapes
:
T

Slice/sizeConst*
valueB:*
dtype0*
_output_shapes
:
b
SliceSliceShape_2Slice/begin
Slice/size*
T0*
Index0*
_output_shapes
:
b
concat/values_0Const*
valueB:
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
:
M
concat/axisConst*
value	B : *
dtype0*
_output_shapes
: 
q
concatConcatV2concat/values_0Sliceconcat/axis*
N*
T0*

Tidx0*
_output_shapes
:
f
	Reshape_2ReshapeLeakyRelu/Maximumconcat*
T0*
Tshape0*
_output_shapes

:2

F
SoftmaxSoftmax	Reshape_2*
T0*
_output_shapes

:2

e
	Reshape_3ReshapeSoftmaxShape_1*
T0*
Tshape0*&
_output_shapes
:2


MatMul_1BatchMatMul	Reshape_3layerfilter0_edgefea_0/Relu*
T0*
adj_x( *
adj_y( *&
_output_shapes
:2 

0BiasAdd/biases/Initializer/zeros/shape_as_tensorConst*
valueB: *
dtype0*!
_class
loc:@BiasAdd/biases*
_output_shapes
:

&BiasAdd/biases/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*!
_class
loc:@BiasAdd/biases*
_output_shapes
: 
Ü
 BiasAdd/biases/Initializer/zerosFill0BiasAdd/biases/Initializer/zeros/shape_as_tensor&BiasAdd/biases/Initializer/zeros/Const*
T0*

index_type0*!
_class
loc:@BiasAdd/biases*
_output_shapes
: 

BiasAdd/biases
VariableV2*
shape: *
dtype0*
	container *
shared_name *!
_class
loc:@BiasAdd/biases*
_output_shapes
: 
Â
BiasAdd/biases/AssignAssignBiasAdd/biases BiasAdd/biases/Initializer/zeros*
T0*
validate_shape(*
use_locking(*!
_class
loc:@BiasAdd/biases*
_output_shapes
: 
w
BiasAdd/biases/readIdentityBiasAdd/biases*
T0*!
_class
loc:@BiasAdd/biases*
_output_shapes
: 

BiasAdd/BiasAddBiasAddMatMul_1BiasAdd/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2 
N
ReluReluBiasAdd/BiasAdd*
T0*&
_output_shapes
:2 
^
concat_1/concat_dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
K
concat_1IdentityRelu*
T0*&
_output_shapes
:2 
[
ExpandDims_6/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
v
ExpandDims_6
ExpandDimsPlaceholderExpandDims_6/dim*
T0*

Tdim0*&
_output_shapes
:2
X
concat_2/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 

concat_2ConcatV2ExpandDims_6concat_1concat_2/axis*
N*
T0*

Tidx0*&
_output_shapes
:2(
^
concat_3/concat_dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
b
concat_3Identitylayerfilter0_edgefea_0/Relu*
T0*&
_output_shapes
:2
 
`
Max/reduction_indicesConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
y
MaxMaxconcat_3Max/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:2 
¯
1gapnet00/weights/Initializer/random_uniform/shapeConst*%
valueB"      (   @   *
dtype0*#
_class
loc:@gapnet00/weights*
_output_shapes
:

/gapnet00/weights/Initializer/random_uniform/minConst*
valueB
 *ôôu¾*
dtype0*#
_class
loc:@gapnet00/weights*
_output_shapes
: 

/gapnet00/weights/Initializer/random_uniform/maxConst*
valueB
 *ôôu>*
dtype0*#
_class
loc:@gapnet00/weights*
_output_shapes
: 
ù
9gapnet00/weights/Initializer/random_uniform/RandomUniformRandomUniform1gapnet00/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@
Þ
/gapnet00/weights/Initializer/random_uniform/subSub/gapnet00/weights/Initializer/random_uniform/max/gapnet00/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet00/weights*
_output_shapes
: 
ø
/gapnet00/weights/Initializer/random_uniform/mulMul9gapnet00/weights/Initializer/random_uniform/RandomUniform/gapnet00/weights/Initializer/random_uniform/sub*
T0*#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@
ê
+gapnet00/weights/Initializer/random_uniformAdd/gapnet00/weights/Initializer/random_uniform/mul/gapnet00/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@
È
gapnet00/weights
VariableV2"/device:CPU:0*
shape:(@*
dtype0*
	container *
shared_name *#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@
î
gapnet00/weights/AssignAssigngapnet00/weights+gapnet00/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@

gapnet00/weights/readIdentitygapnet00/weights"/device:CPU:0*
T0*#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@
Q
gapnet00/L2LossL2Lossgapnet00/weights/read*
T0*
_output_shapes
: 
[
gapnet00/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
e
gapnet00/weight_lossMulgapnet00/L2Lossgapnet00/weight_loss/y*
T0*
_output_shapes
: 
Ú
gapnet00/Conv2DConv2Dconcat_2gapnet00/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2@

!gapnet00/biases/Initializer/ConstConst*
valueB@*    *
dtype0*"
_class
loc:@gapnet00/biases*
_output_shapes
:@
®
gapnet00/biases
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *"
_class
loc:@gapnet00/biases*
_output_shapes
:@
Õ
gapnet00/biases/AssignAssigngapnet00/biases!gapnet00/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet00/biases*
_output_shapes
:@
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
:2@
^
gapnet00/bn/ConstConst*
valueB@*    *
dtype0*
_output_shapes
:@
|
gapnet00/bn/beta
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@
¹
gapnet00/bn/beta/AssignAssigngapnet00/bn/betagapnet00/bn/Const*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet00/bn/beta*
_output_shapes
:@
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
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@
¾
gapnet00/bn/gamma/AssignAssigngapnet00/bn/gammagapnet00/bn/Const_1*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet00/bn/gamma*
_output_shapes
:@
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
gapnet00/bn/moments/meanMeangapnet00/BiasAdd*gapnet00/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@
{
 gapnet00/bn/moments/StopGradientStopGradientgapnet00/bn/moments/mean*
T0*&
_output_shapes
:@

%gapnet00/bn/moments/SquaredDifferenceSquaredDifferencegapnet00/BiasAdd gapnet00/bn/moments/StopGradient*
T0*&
_output_shapes
:2@

.gapnet00/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
É
gapnet00/bn/moments/varianceMean%gapnet00/bn/moments/SquaredDifference.gapnet00/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@
~
gapnet00/bn/moments/SqueezeSqueezegapnet00/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:@

gapnet00/bn/moments/Squeeze_1Squeezegapnet00/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:@
b
gapnet00/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
a
gapnet00/bn/cond/switch_tIdentitygapnet00/bn/cond/Switch:1*
T0
*
_output_shapes
: 
_
gapnet00/bn/cond/switch_fIdentitygapnet00/bn/cond/Switch*
T0
*
_output_shapes
: 
T
gapnet00/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

bgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ò
Xgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¤
Rgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFillbgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorXgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageRgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

Egapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

dgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ö
Zgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¬
Tgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilldgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorZgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Bgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageTgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Ggapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

/gapnet00/bn/cond/ExponentialMovingAverage/decayConst^gapnet00/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
õ
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet00/bn/cond/switch_t*
valueB
 *  ?*
dtype0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¬
=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x/gapnet00/bn/cond/ExponentialMovingAverage/decay*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ö
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubHgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
±
Fgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchEgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet00/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
ä
Hgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet00/bn/moments/Squeezegapnet00/bn/cond/pred_id*
T0*.
_class$
" loc:@gapnet00/bn/moments/Squeeze* 
_output_shapes
:@:@
¾
=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ö
9gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubBgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1=gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *S
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
ù
Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet00/bn/cond/switch_t*
valueB
 *  ?*
dtype0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
²
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubAgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x/gapnet00/bn/cond/ExponentialMovingAverage/decay*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Þ
Agapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubJgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Lgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
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
Æ
?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulAgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Þ
;gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubDgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1?gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
¯
Bgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet00/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Ç
)gapnet00/bn/cond/ExponentialMovingAverageNoOp:^gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg<^gapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1^gapnet00/bn/cond/switch_t
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
%gapnet00/bn/cond/control_dependency_1Identitygapnet00/bn/cond/switch_f^gapnet00/bn/cond/NoOp*
T0
*,
_class"
 loc:@gapnet00/bn/cond/switch_f*
_output_shapes
: 

gapnet00/bn/cond/MergeMerge%gapnet00/bn/cond/control_dependency_1#gapnet00/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
d
gapnet00/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
e
gapnet00/bn/cond_1/switch_tIdentitygapnet00/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
c
gapnet00/bn/cond_1/switch_fIdentitygapnet00/bn/cond_1/Switch*
T0
*
_output_shapes
: 
V
gapnet00/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

gapnet00/bn/cond_1/IdentityIdentity$gapnet00/bn/cond_1/Identity/Switch:1^gapnet00/bn/cond/Merge*
T0*
_output_shapes
:@
À
"gapnet00/bn/cond_1/Identity/SwitchSwitchgapnet00/bn/moments/Squeezegapnet00/bn/cond_1/pred_id*
T0*.
_class$
" loc:@gapnet00/bn/moments/Squeeze* 
_output_shapes
:@:@

gapnet00/bn/cond_1/Identity_1Identity&gapnet00/bn/cond_1/Identity_1/Switch:1^gapnet00/bn/cond/Merge*
T0*
_output_shapes
:@
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
T0*
N*
_output_shapes

:@: 
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
gapnet00/bn/batchnorm/RsqrtRsqrtgapnet00/bn/batchnorm/add*
T0*
_output_shapes
:@
z
gapnet00/bn/batchnorm/mulMulgapnet00/bn/batchnorm/Rsqrtgapnet00/bn/gamma/read*
T0*
_output_shapes
:@

gapnet00/bn/batchnorm/mul_1Mulgapnet00/BiasAddgapnet00/bn/batchnorm/mul*
T0*&
_output_shapes
:2@
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
gapnet00/bn/batchnorm/add_1Addgapnet00/bn/batchnorm/mul_1gapnet00/bn/batchnorm/sub*
T0*&
_output_shapes
:2@
c
gapnet00/ReluRelugapnet00/bn/batchnorm/add_1*
T0*&
_output_shapes
:2@
¯
1gapnet01/weights/Initializer/random_uniform/shapeConst*%
valueB"      @      *
dtype0*#
_class
loc:@gapnet01/weights*
_output_shapes
:

/gapnet01/weights/Initializer/random_uniform/minConst*
valueB
 *ó5¾*
dtype0*#
_class
loc:@gapnet01/weights*
_output_shapes
: 

/gapnet01/weights/Initializer/random_uniform/maxConst*
valueB
 *ó5>*
dtype0*#
_class
loc:@gapnet01/weights*
_output_shapes
: 
ú
9gapnet01/weights/Initializer/random_uniform/RandomUniformRandomUniform1gapnet01/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*#
_class
loc:@gapnet01/weights*'
_output_shapes
:@
Þ
/gapnet01/weights/Initializer/random_uniform/subSub/gapnet01/weights/Initializer/random_uniform/max/gapnet01/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet01/weights*
_output_shapes
: 
ù
/gapnet01/weights/Initializer/random_uniform/mulMul9gapnet01/weights/Initializer/random_uniform/RandomUniform/gapnet01/weights/Initializer/random_uniform/sub*
T0*#
_class
loc:@gapnet01/weights*'
_output_shapes
:@
ë
+gapnet01/weights/Initializer/random_uniformAdd/gapnet01/weights/Initializer/random_uniform/mul/gapnet01/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet01/weights*'
_output_shapes
:@
Ê
gapnet01/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *#
_class
loc:@gapnet01/weights*'
_output_shapes
:@
ï
gapnet01/weights/AssignAssigngapnet01/weights+gapnet01/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet01/weights*'
_output_shapes
:@

gapnet01/weights/readIdentitygapnet01/weights"/device:CPU:0*
T0*#
_class
loc:@gapnet01/weights*'
_output_shapes
:@
Q
gapnet01/L2LossL2Lossgapnet01/weights/read*
T0*
_output_shapes
: 
[
gapnet01/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
e
gapnet01/weight_lossMulgapnet01/L2Lossgapnet01/weight_loss/y*
T0*
_output_shapes
: 
à
gapnet01/Conv2DConv2Dgapnet00/Relugapnet01/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*'
_output_shapes
:2

!gapnet01/biases/Initializer/ConstConst*
valueB*    *
dtype0*"
_class
loc:@gapnet01/biases*
_output_shapes	
:
°
gapnet01/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *"
_class
loc:@gapnet01/biases*
_output_shapes	
:
Ö
gapnet01/biases/AssignAssigngapnet01/biases!gapnet01/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet01/biases*
_output_shapes	
:

gapnet01/biases/readIdentitygapnet01/biases"/device:CPU:0*
T0*"
_class
loc:@gapnet01/biases*
_output_shapes	
:

gapnet01/BiasAddBiasAddgapnet01/Conv2Dgapnet01/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:2
`
gapnet01/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
~
gapnet01/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
º
gapnet01/bn/beta/AssignAssigngapnet01/bn/betagapnet01/bn/Const*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet01/bn/beta*
_output_shapes	
:
~
gapnet01/bn/beta/readIdentitygapnet01/bn/beta*
T0*#
_class
loc:@gapnet01/bn/beta*
_output_shapes	
:
b
gapnet01/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:

gapnet01/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
¿
gapnet01/bn/gamma/AssignAssigngapnet01/bn/gammagapnet01/bn/Const_1*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet01/bn/gamma*
_output_shapes	
:

gapnet01/bn/gamma/readIdentitygapnet01/bn/gamma*
T0*$
_class
loc:@gapnet01/bn/gamma*
_output_shapes	
:

*gapnet01/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
­
gapnet01/bn/moments/meanMeangapnet01/BiasAdd*gapnet01/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:
|
 gapnet01/bn/moments/StopGradientStopGradientgapnet01/bn/moments/mean*
T0*'
_output_shapes
:
 
%gapnet01/bn/moments/SquaredDifferenceSquaredDifferencegapnet01/BiasAdd gapnet01/bn/moments/StopGradient*
T0*'
_output_shapes
:2

.gapnet01/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ê
gapnet01/bn/moments/varianceMean%gapnet01/bn/moments/SquaredDifference.gapnet01/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:

gapnet01/bn/moments/SqueezeSqueezegapnet01/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes	
:

gapnet01/bn/moments/Squeeze_1Squeezegapnet01/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes	
:
b
gapnet01/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
a
gapnet01/bn/cond/switch_tIdentitygapnet01/bn/cond/Switch:1*
T0
*
_output_shapes
: 
_
gapnet01/bn/cond/switch_fIdentitygapnet01/bn/cond/Switch*
T0
*
_output_shapes
: 
T
gapnet01/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

bgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ò
Xgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¥
Rgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFillbgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorXgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Ggapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverageRgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Egapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
T0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

dgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ö
Zgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
­
Tgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilldgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorZgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Bgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Igapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverageTgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Ggapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

/gapnet01/bn/cond/ExponentialMovingAverage/decayConst^gapnet01/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
õ
?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet01/bn/cond/switch_t*
valueB
 *  ?*
dtype0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¬
=gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x/gapnet01/bn/cond/ExponentialMovingAverage/decay*
T0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
×
?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubHgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Jgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
³
Fgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchEgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet01/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
æ
Hgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet01/bn/moments/Squeezegapnet01/bn/cond/pred_id*
T0*.
_class$
" loc:@gapnet01/bn/moments/Squeeze*"
_output_shapes
::
¿
=gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1=gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
×
9gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubBgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1=gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
«
@gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAveragegapnet01/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
ù
Agapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet01/bn/cond/switch_t*
valueB
 *  ?*
dtype0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
²
?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubAgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x/gapnet01/bn/cond/ExponentialMovingAverage/decay*
T0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ß
Agapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubJgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Lgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
¹
Hgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchGgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet01/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
ì
Jgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchgapnet01/bn/moments/Squeeze_1gapnet01/bn/cond/pred_id*
T0*0
_class&
$"loc:@gapnet01/bn/moments/Squeeze_1*"
_output_shapes
::
Ç
?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulAgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
ß
;gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubDgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
±
Bgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet01/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
Ç
)gapnet01/bn/cond/ExponentialMovingAverageNoOp:^gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg<^gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1^gapnet01/bn/cond/switch_t
Å
#gapnet01/bn/cond/control_dependencyIdentitygapnet01/bn/cond/switch_t*^gapnet01/bn/cond/ExponentialMovingAverage*
T0
*,
_class"
 loc:@gapnet01/bn/cond/switch_t*
_output_shapes
: 
9
gapnet01/bn/cond/NoOpNoOp^gapnet01/bn/cond/switch_f
³
%gapnet01/bn/cond/control_dependency_1Identitygapnet01/bn/cond/switch_f^gapnet01/bn/cond/NoOp*
T0
*,
_class"
 loc:@gapnet01/bn/cond/switch_f*
_output_shapes
: 

gapnet01/bn/cond/MergeMerge%gapnet01/bn/cond/control_dependency_1#gapnet01/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
d
gapnet01/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
e
gapnet01/bn/cond_1/switch_tIdentitygapnet01/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
c
gapnet01/bn/cond_1/switch_fIdentitygapnet01/bn/cond_1/Switch*
T0
*
_output_shapes
: 
V
gapnet01/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

gapnet01/bn/cond_1/IdentityIdentity$gapnet01/bn/cond_1/Identity/Switch:1^gapnet01/bn/cond/Merge*
T0*
_output_shapes	
:
Â
"gapnet01/bn/cond_1/Identity/SwitchSwitchgapnet01/bn/moments/Squeezegapnet01/bn/cond_1/pred_id*
T0*.
_class$
" loc:@gapnet01/bn/moments/Squeeze*"
_output_shapes
::

gapnet01/bn/cond_1/Identity_1Identity&gapnet01/bn/cond_1/Identity_1/Switch:1^gapnet01/bn/cond/Merge*
T0*
_output_shapes	
:
È
$gapnet01/bn/cond_1/Identity_1/SwitchSwitchgapnet01/bn/moments/Squeeze_1gapnet01/bn/cond_1/pred_id*
T0*0
_class&
$"loc:@gapnet01/bn/moments/Squeeze_1*"
_output_shapes
::

gapnet01/bn/cond_1/Switch_1SwitchEgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet01/bn/cond_1/pred_id*
T0*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

gapnet01/bn/cond_1/Switch_2SwitchGgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet01/bn/cond_1/pred_id*
T0*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::

gapnet01/bn/cond_1/MergeMergegapnet01/bn/cond_1/Switch_1gapnet01/bn/cond_1/Identity*
T0*
N*
_output_shapes
	:: 

gapnet01/bn/cond_1/Merge_1Mergegapnet01/bn/cond_1/Switch_2gapnet01/bn/cond_1/Identity_1*
T0*
N*
_output_shapes
	:: 
`
gapnet01/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 

gapnet01/bn/batchnorm/addAddgapnet01/bn/cond_1/Merge_1gapnet01/bn/batchnorm/add/y*
T0*
_output_shapes	
:
e
gapnet01/bn/batchnorm/RsqrtRsqrtgapnet01/bn/batchnorm/add*
T0*
_output_shapes	
:
{
gapnet01/bn/batchnorm/mulMulgapnet01/bn/batchnorm/Rsqrtgapnet01/bn/gamma/read*
T0*
_output_shapes	
:

gapnet01/bn/batchnorm/mul_1Mulgapnet01/BiasAddgapnet01/bn/batchnorm/mul*
T0*'
_output_shapes
:2
}
gapnet01/bn/batchnorm/mul_2Mulgapnet01/bn/cond_1/Mergegapnet01/bn/batchnorm/mul*
T0*
_output_shapes	
:
z
gapnet01/bn/batchnorm/subSubgapnet01/bn/beta/readgapnet01/bn/batchnorm/mul_2*
T0*
_output_shapes	
:

gapnet01/bn/batchnorm/add_1Addgapnet01/bn/batchnorm/mul_1gapnet01/bn/batchnorm/sub*
T0*'
_output_shapes
:2
d
gapnet01/ReluRelugapnet01/bn/batchnorm/add_1*
T0*'
_output_shapes
:2
a
	Squeeze_3Squeezegapnet01/Relu*
T0*
squeeze_dims
 *
_output_shapes
:	2
R
ExpandDims_7/dimConst*
value	B : *
dtype0*
_output_shapes
: 
q
ExpandDims_7
ExpandDims	Squeeze_3ExpandDims_7/dim*
T0*

Tdim0*#
_output_shapes
:2
j
strided_slice_5/stackConst*!
valueB"           *
dtype0*
_output_shapes
:
l
strided_slice_5/stack_1Const*!
valueB"           *
dtype0*
_output_shapes
:
l
strided_slice_5/stack_2Const*!
valueB"         *
dtype0*
_output_shapes
:

strided_slice_5StridedSliceExpandDims_7strided_slice_5/stackstrided_slice_5/stack_1strided_slice_5/stack_2*
T0*
Index0*

begin_mask*
end_mask*
ellipsis_mask *
new_axis_mask *
shrink_axis_mask*
_output_shapes

:2
[
ExpandDims_8/dimConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
v
ExpandDims_8
ExpandDimsstrided_slice_5ExpandDims_8/dim*
T0*

Tdim0*"
_output_shapes
:2
e
transpose_5/permConst*!
valueB"          *
dtype0*
_output_shapes
:
s
transpose_5	TransposeExpandDims_7transpose_5/perm*
T0*
Tperm0*#
_output_shapes
:2
y
MatMul_2BatchMatMulExpandDims_7transpose_5*
T0*
adj_x( *
adj_y( *"
_output_shapes
:22
L
mul_4/xConst*
valueB
 *   À*
dtype0*
_output_shapes
: 
L
mul_4Mulmul_4/xMatMul_2*
T0*"
_output_shapes
:22
N
Square_1SquareExpandDims_7*
T0*#
_output_shapes
:2
b
Sum_1/reduction_indicesConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
y
Sum_1SumSquare_1Sum_1/reduction_indices*
	keep_dims(*
T0*

Tidx0*"
_output_shapes
:2
e
transpose_6/permConst*!
valueB"          *
dtype0*
_output_shapes
:
k
transpose_6	TransposeSum_1transpose_6/perm*
T0*
Tperm0*"
_output_shapes
:2
G
add_7AddSum_1mul_4*
T0*"
_output_shapes
:22
M
add_8Addadd_7transpose_6*
T0*"
_output_shapes
:22
@
Neg_1Negadd_8*
T0*"
_output_shapes
:22
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
:2
:2

a
	Squeeze_4Squeezegapnet01/Relu*
T0*
squeeze_dims
 *
_output_shapes
:	2
R
ExpandDims_9/dimConst*
value	B : *
dtype0*
_output_shapes
: 
q
ExpandDims_9
ExpandDims	Squeeze_4ExpandDims_9/dim*
T0*

Tdim0*#
_output_shapes
:2
\
ExpandDims_10/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
z
ExpandDims_10
ExpandDimsExpandDims_9ExpandDims_10/dim*
T0*

Tdim0*'
_output_shapes
:2
Ý
Hlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"         @   *
dtype0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*
_output_shapes
:
Ç
Flayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *ó5¾*
dtype0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*
_output_shapes
: 
Ç
Flayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *ó5>*
dtype0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*
_output_shapes
: 
¿
Playerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformHlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@
º
Flayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/subSubFlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/maxFlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/min*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*
_output_shapes
: 
Õ
Flayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/mulMulPlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/RandomUniformFlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/sub*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@
Ç
Blayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniformAddFlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/mulFlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform/min*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@
ø
'layerfilter1_newfea_conv_head_0/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@
Ë
.layerfilter1_newfea_conv_head_0/weights/AssignAssign'layerfilter1_newfea_conv_head_0/weightsBlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@
Þ
,layerfilter1_newfea_conv_head_0/weights/readIdentity'layerfilter1_newfea_conv_head_0/weights"/device:CPU:0*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@

&layerfilter1_newfea_conv_head_0/L2LossL2Loss,layerfilter1_newfea_conv_head_0/weights/read*
T0*
_output_shapes
: 
r
-layerfilter1_newfea_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
ª
+layerfilter1_newfea_conv_head_0/weight_lossMul&layerfilter1_newfea_conv_head_0/L2Loss-layerfilter1_newfea_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 

&layerfilter1_newfea_conv_head_0/Conv2DConv2DExpandDims_10,layerfilter1_newfea_conv_head_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2@
u
(layerfilter1_newfea_conv_head_0/bn/ConstConst*
valueB@*    *
dtype0*
_output_shapes
:@

'layerfilter1_newfea_conv_head_0/bn/beta
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@

.layerfilter1_newfea_conv_head_0/bn/beta/AssignAssign'layerfilter1_newfea_conv_head_0/bn/beta(layerfilter1_newfea_conv_head_0/bn/Const*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/bn/beta*
_output_shapes
:@
Â
,layerfilter1_newfea_conv_head_0/bn/beta/readIdentity'layerfilter1_newfea_conv_head_0/bn/beta*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/bn/beta*
_output_shapes
:@
w
*layerfilter1_newfea_conv_head_0/bn/Const_1Const*
valueB@*  ?*
dtype0*
_output_shapes
:@

(layerfilter1_newfea_conv_head_0/bn/gamma
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@

/layerfilter1_newfea_conv_head_0/bn/gamma/AssignAssign(layerfilter1_newfea_conv_head_0/bn/gamma*layerfilter1_newfea_conv_head_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_newfea_conv_head_0/bn/gamma*
_output_shapes
:@
Å
-layerfilter1_newfea_conv_head_0/bn/gamma/readIdentity(layerfilter1_newfea_conv_head_0/bn/gamma*
T0*;
_class1
/-loc:@layerfilter1_newfea_conv_head_0/bn/gamma*
_output_shapes
:@

Alayerfilter1_newfea_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ð
/layerfilter1_newfea_conv_head_0/bn/moments/meanMean&layerfilter1_newfea_conv_head_0/Conv2DAlayerfilter1_newfea_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@
©
7layerfilter1_newfea_conv_head_0/bn/moments/StopGradientStopGradient/layerfilter1_newfea_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
:@
ã
<layerfilter1_newfea_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference&layerfilter1_newfea_conv_head_0/Conv2D7layerfilter1_newfea_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2@

Elayerfilter1_newfea_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

3layerfilter1_newfea_conv_head_0/bn/moments/varianceMean<layerfilter1_newfea_conv_head_0/bn/moments/SquaredDifferenceElayerfilter1_newfea_conv_head_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@
¬
2layerfilter1_newfea_conv_head_0/bn/moments/SqueezeSqueeze/layerfilter1_newfea_conv_head_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:@
²
4layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1Squeeze3layerfilter1_newfea_conv_head_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:@
y
.layerfilter1_newfea_conv_head_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

0layerfilter1_newfea_conv_head_0/bn/cond/switch_tIdentity0layerfilter1_newfea_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

0layerfilter1_newfea_conv_head_0/bn/cond/switch_fIdentity.layerfilter1_newfea_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
k
/layerfilter1_newfea_conv_head_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ß
layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ð
layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
à
layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Þ
nlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ä
ulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

slayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentitynlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ã
layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ô
layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
è
layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
â
playerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Ì
wlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

ulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
¾
Flayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst1^layerfilter1_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ò
Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst1^layerfilter1_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
 
Tlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubVlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xFlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ê
Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Sub_layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1alayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
¼
]layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchslayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read/layerfilter1_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
À
_layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch2layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/layerfilter1_newfea_conv_head_0/bn/cond/pred_id*
T0*E
_class;
97loc:@layerfilter1_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
:@:@
²
Tlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulVlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Tlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ê
Playerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubYlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Tlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
´
Wlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/layerfilter1_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
Ö
Xlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst1^layerfilter1_newfea_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¦
Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubXlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xFlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ò
Xlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subalayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1clayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Â
_layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read/layerfilter1_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Æ
alayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch4layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/layerfilter1_newfea_conv_head_0/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
:@:@
º
Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulXlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Ò
Rlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub[layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
º
Ylayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/layerfilter1_newfea_conv_head_0/bn/cond/pred_id*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
£
@layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverageNoOpQ^layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgS^layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_11^layerfilter1_newfea_conv_head_0/bn/cond/switch_t
¡
:layerfilter1_newfea_conv_head_0/bn/cond/control_dependencyIdentity0layerfilter1_newfea_conv_head_0/bn/cond/switch_tA^layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*C
_class9
75loc:@layerfilter1_newfea_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
g
,layerfilter1_newfea_conv_head_0/bn/cond/NoOpNoOp1^layerfilter1_newfea_conv_head_0/bn/cond/switch_f

<layerfilter1_newfea_conv_head_0/bn/cond/control_dependency_1Identity0layerfilter1_newfea_conv_head_0/bn/cond/switch_f-^layerfilter1_newfea_conv_head_0/bn/cond/NoOp*
T0
*C
_class9
75loc:@layerfilter1_newfea_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
Ü
-layerfilter1_newfea_conv_head_0/bn/cond/MergeMerge<layerfilter1_newfea_conv_head_0/bn/cond/control_dependency_1:layerfilter1_newfea_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
{
0layerfilter1_newfea_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter1_newfea_conv_head_0/bn/cond_1/switch_tIdentity2layerfilter1_newfea_conv_head_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

2layerfilter1_newfea_conv_head_0/bn/cond_1/switch_fIdentity0layerfilter1_newfea_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
m
1layerfilter1_newfea_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ð
2layerfilter1_newfea_conv_head_0/bn/cond_1/IdentityIdentity;layerfilter1_newfea_conv_head_0/bn/cond_1/Identity/Switch:1.^layerfilter1_newfea_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:@

9layerfilter1_newfea_conv_head_0/bn/cond_1/Identity/SwitchSwitch2layerfilter1_newfea_conv_head_0/bn/moments/Squeeze1layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id*
T0*E
_class;
97loc:@layerfilter1_newfea_conv_head_0/bn/moments/Squeeze* 
_output_shapes
:@:@
Ô
4layerfilter1_newfea_conv_head_0/bn/cond_1/Identity_1Identity=layerfilter1_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1.^layerfilter1_newfea_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:@
¢
;layerfilter1_newfea_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch4layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_11layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
:@:@

2layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_1Switchslayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@

2layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_2Switchulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Ð
/layerfilter1_newfea_conv_head_0/bn/cond_1/MergeMerge2layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_12layerfilter1_newfea_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:@: 
Ô
1layerfilter1_newfea_conv_head_0/bn/cond_1/Merge_1Merge2layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_24layerfilter1_newfea_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:@: 
w
2layerfilter1_newfea_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
Ã
0layerfilter1_newfea_conv_head_0/bn/batchnorm/addAdd1layerfilter1_newfea_conv_head_0/bn/cond_1/Merge_12layerfilter1_newfea_conv_head_0/bn/batchnorm/add/y*
T0*
_output_shapes
:@

2layerfilter1_newfea_conv_head_0/bn/batchnorm/RsqrtRsqrt0layerfilter1_newfea_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
:@
¿
0layerfilter1_newfea_conv_head_0/bn/batchnorm/mulMul2layerfilter1_newfea_conv_head_0/bn/batchnorm/Rsqrt-layerfilter1_newfea_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:@
Ä
2layerfilter1_newfea_conv_head_0/bn/batchnorm/mul_1Mul&layerfilter1_newfea_conv_head_0/Conv2D0layerfilter1_newfea_conv_head_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2@
Á
2layerfilter1_newfea_conv_head_0/bn/batchnorm/mul_2Mul/layerfilter1_newfea_conv_head_0/bn/cond_1/Merge0layerfilter1_newfea_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:@
¾
0layerfilter1_newfea_conv_head_0/bn/batchnorm/subSub,layerfilter1_newfea_conv_head_0/bn/beta/read2layerfilter1_newfea_conv_head_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:@
Ð
2layerfilter1_newfea_conv_head_0/bn/batchnorm/add_1Add2layerfilter1_newfea_conv_head_0/bn/batchnorm/mul_10layerfilter1_newfea_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2@

$layerfilter1_newfea_conv_head_0/ReluRelu2layerfilter1_newfea_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2@
a
	Squeeze_5SqueezeExpandDims_10*
T0*
squeeze_dims
 *
_output_shapes
:	2
S
ExpandDims_11/dimConst*
value	B : *
dtype0*
_output_shapes
: 
s
ExpandDims_11
ExpandDims	Squeeze_5ExpandDims_11/dim*
T0*

Tdim0*#
_output_shapes
:2
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
value	B :2*
dtype0*
_output_shapes
: 
C
mul_5Mulrange_1mul_5/y*
T0*
_output_shapes
:
d
Reshape_4/shapeConst*!
valueB"         *
dtype0*
_output_shapes
:
g
	Reshape_4Reshapemul_5Reshape_4/shape*
T0*
Tshape0*"
_output_shapes
:
`
Reshape_5/shapeConst*
valueB"ÿÿÿÿ   *
dtype0*
_output_shapes
:
l
	Reshape_5ReshapeExpandDims_11Reshape_5/shape*
T0*
Tshape0*
_output_shapes
:	2
P
add_9Add
TopKV2_1:1	Reshape_4*
T0*"
_output_shapes
:2


Gather_1Gather	Reshape_5add_9*
validate_indices(*
Tparams0*
Tindices0*'
_output_shapes
:2

i
Tile_2/multiplesConst*%
valueB"      
      *
dtype0*
_output_shapes
:
s
Tile_2TileExpandDims_10Tile_2/multiples*
T0*

Tmultiples0*'
_output_shapes
:2

P
sub_4SubTile_2Gather_1*
T0*'
_output_shapes
:2

Ë
?layerfilter1_edgefea_0/weights/Initializer/random_uniform/shapeConst*%
valueB"         @   *
dtype0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*
_output_shapes
:
µ
=layerfilter1_edgefea_0/weights/Initializer/random_uniform/minConst*
valueB
 *ó5¾*
dtype0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*
_output_shapes
: 
µ
=layerfilter1_edgefea_0/weights/Initializer/random_uniform/maxConst*
valueB
 *ó5>*
dtype0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*
_output_shapes
: 
¤
Glayerfilter1_edgefea_0/weights/Initializer/random_uniform/RandomUniformRandomUniform?layerfilter1_edgefea_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@

=layerfilter1_edgefea_0/weights/Initializer/random_uniform/subSub=layerfilter1_edgefea_0/weights/Initializer/random_uniform/max=layerfilter1_edgefea_0/weights/Initializer/random_uniform/min*
T0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*
_output_shapes
: 
±
=layerfilter1_edgefea_0/weights/Initializer/random_uniform/mulMulGlayerfilter1_edgefea_0/weights/Initializer/random_uniform/RandomUniform=layerfilter1_edgefea_0/weights/Initializer/random_uniform/sub*
T0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@
£
9layerfilter1_edgefea_0/weights/Initializer/random_uniformAdd=layerfilter1_edgefea_0/weights/Initializer/random_uniform/mul=layerfilter1_edgefea_0/weights/Initializer/random_uniform/min*
T0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@
æ
layerfilter1_edgefea_0/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@
§
%layerfilter1_edgefea_0/weights/AssignAssignlayerfilter1_edgefea_0/weights9layerfilter1_edgefea_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@
Ã
#layerfilter1_edgefea_0/weights/readIdentitylayerfilter1_edgefea_0/weights"/device:CPU:0*
T0*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@
m
layerfilter1_edgefea_0/L2LossL2Loss#layerfilter1_edgefea_0/weights/read*
T0*
_output_shapes
: 
i
$layerfilter1_edgefea_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 

"layerfilter1_edgefea_0/weight_lossMullayerfilter1_edgefea_0/L2Loss$layerfilter1_edgefea_0/weight_loss/y*
T0*
_output_shapes
: 
ó
layerfilter1_edgefea_0/Conv2DConv2Dsub_4#layerfilter1_edgefea_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2
@
®
/layerfilter1_edgefea_0/biases/Initializer/ConstConst*
valueB@*    *
dtype0*0
_class&
$"loc:@layerfilter1_edgefea_0/biases*
_output_shapes
:@
Ê
layerfilter1_edgefea_0/biases
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *0
_class&
$"loc:@layerfilter1_edgefea_0/biases*
_output_shapes
:@

$layerfilter1_edgefea_0/biases/AssignAssignlayerfilter1_edgefea_0/biases/layerfilter1_edgefea_0/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter1_edgefea_0/biases*
_output_shapes
:@
³
"layerfilter1_edgefea_0/biases/readIdentitylayerfilter1_edgefea_0/biases"/device:CPU:0*
T0*0
_class&
$"loc:@layerfilter1_edgefea_0/biases*
_output_shapes
:@
´
layerfilter1_edgefea_0/BiasAddBiasAddlayerfilter1_edgefea_0/Conv2D"layerfilter1_edgefea_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
@
l
layerfilter1_edgefea_0/bn/ConstConst*
valueB@*    *
dtype0*
_output_shapes
:@

layerfilter1_edgefea_0/bn/beta
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@
ñ
%layerfilter1_edgefea_0/bn/beta/AssignAssignlayerfilter1_edgefea_0/bn/betalayerfilter1_edgefea_0/bn/Const*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_0/bn/beta*
_output_shapes
:@
§
#layerfilter1_edgefea_0/bn/beta/readIdentitylayerfilter1_edgefea_0/bn/beta*
T0*1
_class'
%#loc:@layerfilter1_edgefea_0/bn/beta*
_output_shapes
:@
n
!layerfilter1_edgefea_0/bn/Const_1Const*
valueB@*  ?*
dtype0*
_output_shapes
:@

layerfilter1_edgefea_0/bn/gamma
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@
ö
&layerfilter1_edgefea_0/bn/gamma/AssignAssignlayerfilter1_edgefea_0/bn/gamma!layerfilter1_edgefea_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter1_edgefea_0/bn/gamma*
_output_shapes
:@
ª
$layerfilter1_edgefea_0/bn/gamma/readIdentitylayerfilter1_edgefea_0/bn/gamma*
T0*2
_class(
&$loc:@layerfilter1_edgefea_0/bn/gamma*
_output_shapes
:@

8layerfilter1_edgefea_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ö
&layerfilter1_edgefea_0/bn/moments/meanMeanlayerfilter1_edgefea_0/BiasAdd8layerfilter1_edgefea_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@

.layerfilter1_edgefea_0/bn/moments/StopGradientStopGradient&layerfilter1_edgefea_0/bn/moments/mean*
T0*&
_output_shapes
:@
É
3layerfilter1_edgefea_0/bn/moments/SquaredDifferenceSquaredDifferencelayerfilter1_edgefea_0/BiasAdd.layerfilter1_edgefea_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2
@

<layerfilter1_edgefea_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ó
*layerfilter1_edgefea_0/bn/moments/varianceMean3layerfilter1_edgefea_0/bn/moments/SquaredDifference<layerfilter1_edgefea_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@

)layerfilter1_edgefea_0/bn/moments/SqueezeSqueeze&layerfilter1_edgefea_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:@
 
+layerfilter1_edgefea_0/bn/moments/Squeeze_1Squeeze*layerfilter1_edgefea_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:@
p
%layerfilter1_edgefea_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
}
'layerfilter1_edgefea_0/bn/cond/switch_tIdentity'layerfilter1_edgefea_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 
{
'layerfilter1_edgefea_0/bn/cond/switch_fIdentity%layerfilter1_edgefea_0/bn/cond/Switch*
T0
*
_output_shapes
: 
b
&layerfilter1_edgefea_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
¹
~layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ª
tlayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

nlayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFill~layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensortlayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
¹
\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ú
clayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragenlayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
á
alayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
¾
layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
®
vlayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

playerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorvlayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
½
^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

elayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssign^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageplayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
ç
clayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentity^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
¬
=layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/decayConst(^layerfilter1_edgefea_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
­
Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst(^layerfilter1_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ò
Klayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubMlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x=layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
T0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubVlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Xlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

Tlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchalayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read&layerfilter1_edgefea_0/bn/cond/pred_id*
T0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@

Vlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch)layerfilter1_edgefea_0/bn/moments/Squeeze&layerfilter1_edgefea_0/bn/cond/pred_id*
T0*<
_class2
0.loc:@layerfilter1_edgefea_0/bn/moments/Squeeze* 
_output_shapes
:@:@

Klayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulMlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Klayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

Glayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubPlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Klayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ý
Nlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage&layerfilter1_edgefea_0/bn/cond/pred_id*
T0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
±
Olayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst(^layerfilter1_edgefea_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubOlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x=layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/decay*
T0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¤
Olayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubXlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Zlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Vlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchclayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read&layerfilter1_edgefea_0/bn/cond/pred_id*
T0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
¢
Xlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch+layerfilter1_edgefea_0/bn/moments/Squeeze_1&layerfilter1_edgefea_0/bn/cond/pred_id*
T0*>
_class4
20loc:@layerfilter1_edgefea_0/bn/moments/Squeeze_1* 
_output_shapes
:@:@

Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulOlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
¤
Ilayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubRlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Playerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitch^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage&layerfilter1_edgefea_0/bn/cond/pred_id*
T0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
ÿ
7layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverageNoOpH^layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgJ^layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1(^layerfilter1_edgefea_0/bn/cond/switch_t
ý
1layerfilter1_edgefea_0/bn/cond/control_dependencyIdentity'layerfilter1_edgefea_0/bn/cond/switch_t8^layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage*
T0
*:
_class0
.,loc:@layerfilter1_edgefea_0/bn/cond/switch_t*
_output_shapes
: 
U
#layerfilter1_edgefea_0/bn/cond/NoOpNoOp(^layerfilter1_edgefea_0/bn/cond/switch_f
ë
3layerfilter1_edgefea_0/bn/cond/control_dependency_1Identity'layerfilter1_edgefea_0/bn/cond/switch_f$^layerfilter1_edgefea_0/bn/cond/NoOp*
T0
*:
_class0
.,loc:@layerfilter1_edgefea_0/bn/cond/switch_f*
_output_shapes
: 
Á
$layerfilter1_edgefea_0/bn/cond/MergeMerge3layerfilter1_edgefea_0/bn/cond/control_dependency_11layerfilter1_edgefea_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
r
'layerfilter1_edgefea_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

)layerfilter1_edgefea_0/bn/cond_1/switch_tIdentity)layerfilter1_edgefea_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

)layerfilter1_edgefea_0/bn/cond_1/switch_fIdentity'layerfilter1_edgefea_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
d
(layerfilter1_edgefea_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
µ
)layerfilter1_edgefea_0/bn/cond_1/IdentityIdentity2layerfilter1_edgefea_0/bn/cond_1/Identity/Switch:1%^layerfilter1_edgefea_0/bn/cond/Merge*
T0*
_output_shapes
:@
ø
0layerfilter1_edgefea_0/bn/cond_1/Identity/SwitchSwitch)layerfilter1_edgefea_0/bn/moments/Squeeze(layerfilter1_edgefea_0/bn/cond_1/pred_id*
T0*<
_class2
0.loc:@layerfilter1_edgefea_0/bn/moments/Squeeze* 
_output_shapes
:@:@
¹
+layerfilter1_edgefea_0/bn/cond_1/Identity_1Identity4layerfilter1_edgefea_0/bn/cond_1/Identity_1/Switch:1%^layerfilter1_edgefea_0/bn/cond/Merge*
T0*
_output_shapes
:@
þ
2layerfilter1_edgefea_0/bn/cond_1/Identity_1/SwitchSwitch+layerfilter1_edgefea_0/bn/moments/Squeeze_1(layerfilter1_edgefea_0/bn/cond_1/pred_id*
T0*>
_class4
20loc:@layerfilter1_edgefea_0/bn/moments/Squeeze_1* 
_output_shapes
:@:@
Ü
)layerfilter1_edgefea_0/bn/cond_1/Switch_1Switchalayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read(layerfilter1_edgefea_0/bn/cond_1/pred_id*
T0*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
à
)layerfilter1_edgefea_0/bn/cond_1/Switch_2Switchclayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read(layerfilter1_edgefea_0/bn/cond_1/pred_id*
T0*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
µ
&layerfilter1_edgefea_0/bn/cond_1/MergeMerge)layerfilter1_edgefea_0/bn/cond_1/Switch_1)layerfilter1_edgefea_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:@: 
¹
(layerfilter1_edgefea_0/bn/cond_1/Merge_1Merge)layerfilter1_edgefea_0/bn/cond_1/Switch_2+layerfilter1_edgefea_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:@: 
n
)layerfilter1_edgefea_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
¨
'layerfilter1_edgefea_0/bn/batchnorm/addAdd(layerfilter1_edgefea_0/bn/cond_1/Merge_1)layerfilter1_edgefea_0/bn/batchnorm/add/y*
T0*
_output_shapes
:@

)layerfilter1_edgefea_0/bn/batchnorm/RsqrtRsqrt'layerfilter1_edgefea_0/bn/batchnorm/add*
T0*
_output_shapes
:@
¤
'layerfilter1_edgefea_0/bn/batchnorm/mulMul)layerfilter1_edgefea_0/bn/batchnorm/Rsqrt$layerfilter1_edgefea_0/bn/gamma/read*
T0*
_output_shapes
:@
ª
)layerfilter1_edgefea_0/bn/batchnorm/mul_1Mullayerfilter1_edgefea_0/BiasAdd'layerfilter1_edgefea_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2
@
¦
)layerfilter1_edgefea_0/bn/batchnorm/mul_2Mul&layerfilter1_edgefea_0/bn/cond_1/Merge'layerfilter1_edgefea_0/bn/batchnorm/mul*
T0*
_output_shapes
:@
£
'layerfilter1_edgefea_0/bn/batchnorm/subSub#layerfilter1_edgefea_0/bn/beta/read)layerfilter1_edgefea_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:@
µ
)layerfilter1_edgefea_0/bn/batchnorm/add_1Add)layerfilter1_edgefea_0/bn/batchnorm/mul_1'layerfilter1_edgefea_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2
@

layerfilter1_edgefea_0/ReluRelu)layerfilter1_edgefea_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2
@
á
Jlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"      @      *
dtype0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*
_output_shapes
:
Ë
Hlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *¾*
dtype0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*
_output_shapes
: 
Ë
Hlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *>*
dtype0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*
_output_shapes
: 
Ä
Rlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@
Â
Hlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/subSubHlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/maxHlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*
_output_shapes
: 
Ü
Hlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/mulMulRlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformHlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/sub*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@
Î
Dlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniformAddHlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/mulHlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@
ú
)layerfilter1_self_att_conv_head_0/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@
Ò
0layerfilter1_self_att_conv_head_0/weights/AssignAssign)layerfilter1_self_att_conv_head_0/weightsDlayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@
ã
.layerfilter1_self_att_conv_head_0/weights/readIdentity)layerfilter1_self_att_conv_head_0/weights"/device:CPU:0*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@

(layerfilter1_self_att_conv_head_0/L2LossL2Loss.layerfilter1_self_att_conv_head_0/weights/read*
T0*
_output_shapes
: 
t
/layerfilter1_self_att_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter1_self_att_conv_head_0/weight_lossMul(layerfilter1_self_att_conv_head_0/L2Loss/layerfilter1_self_att_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 
¨
(layerfilter1_self_att_conv_head_0/Conv2DConv2D$layerfilter1_newfea_conv_head_0/Relu.layerfilter1_self_att_conv_head_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2
Ä
:layerfilter1_self_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *
dtype0*;
_class1
/-loc:@layerfilter1_self_att_conv_head_0/biases*
_output_shapes
:
à
(layerfilter1_self_att_conv_head_0/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *;
_class1
/-loc:@layerfilter1_self_att_conv_head_0/biases*
_output_shapes
:
¹
/layerfilter1_self_att_conv_head_0/biases/AssignAssign(layerfilter1_self_att_conv_head_0/biases:layerfilter1_self_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_self_att_conv_head_0/biases*
_output_shapes
:
Ô
-layerfilter1_self_att_conv_head_0/biases/readIdentity(layerfilter1_self_att_conv_head_0/biases"/device:CPU:0*
T0*;
_class1
/-loc:@layerfilter1_self_att_conv_head_0/biases*
_output_shapes
:
Õ
)layerfilter1_self_att_conv_head_0/BiasAddBiasAdd(layerfilter1_self_att_conv_head_0/Conv2D-layerfilter1_self_att_conv_head_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
w
*layerfilter1_self_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

)layerfilter1_self_att_conv_head_0/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:

0layerfilter1_self_att_conv_head_0/bn/beta/AssignAssign)layerfilter1_self_att_conv_head_0/bn/beta*layerfilter1_self_att_conv_head_0/bn/Const*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/bn/beta*
_output_shapes
:
È
.layerfilter1_self_att_conv_head_0/bn/beta/readIdentity)layerfilter1_self_att_conv_head_0/bn/beta*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/bn/beta*
_output_shapes
:
y
,layerfilter1_self_att_conv_head_0/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

*layerfilter1_self_att_conv_head_0/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
¢
1layerfilter1_self_att_conv_head_0/bn/gamma/AssignAssign*layerfilter1_self_att_conv_head_0/bn/gamma,layerfilter1_self_att_conv_head_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_self_att_conv_head_0/bn/gamma*
_output_shapes
:
Ë
/layerfilter1_self_att_conv_head_0/bn/gamma/readIdentity*layerfilter1_self_att_conv_head_0/bn/gamma*
T0*=
_class3
1/loc:@layerfilter1_self_att_conv_head_0/bn/gamma*
_output_shapes
:

Clayerfilter1_self_att_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
÷
1layerfilter1_self_att_conv_head_0/bn/moments/meanMean)layerfilter1_self_att_conv_head_0/BiasAddClayerfilter1_self_att_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
­
9layerfilter1_self_att_conv_head_0/bn/moments/StopGradientStopGradient1layerfilter1_self_att_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
:
ê
>layerfilter1_self_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference)layerfilter1_self_att_conv_head_0/BiasAdd9layerfilter1_self_att_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2

Glayerfilter1_self_att_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

5layerfilter1_self_att_conv_head_0/bn/moments/varianceMean>layerfilter1_self_att_conv_head_0/bn/moments/SquaredDifferenceGlayerfilter1_self_att_conv_head_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
°
4layerfilter1_self_att_conv_head_0/bn/moments/SqueezeSqueeze1layerfilter1_self_att_conv_head_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:
¶
6layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1Squeeze5layerfilter1_self_att_conv_head_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:
{
0layerfilter1_self_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter1_self_att_conv_head_0/bn/cond/switch_tIdentity2layerfilter1_self_att_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter1_self_att_conv_head_0/bn/cond/switch_fIdentity0layerfilter1_self_att_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter1_self_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ç
layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ø
layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ð
layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
æ
rlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
ylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¤
wlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ë
layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ê
tlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
{layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ª
ylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Â
Hlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter1_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ú
Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter1_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ª
Vlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ô
Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
È
_layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter1_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
È
alayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter1_self_att_conv_head_0/bn/moments/Squeeze1layerfilter1_self_att_conv_head_0/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter1_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
¼
Vlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
Rlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
À
Ylayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter1_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Þ
Zlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter1_self_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
°
Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ü
Zlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Î
alayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter1_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Î
clayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_11layerfilter1_self_att_conv_head_0/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
Ä
Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
Tlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
[layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter1_self_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
«
Blayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpS^layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_13^layerfilter1_self_att_conv_head_0/bn/cond/switch_t
©
<layerfilter1_self_att_conv_head_0/bn/cond/control_dependencyIdentity2layerfilter1_self_att_conv_head_0/bn/cond/switch_tC^layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter1_self_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter1_self_att_conv_head_0/bn/cond/NoOpNoOp3^layerfilter1_self_att_conv_head_0/bn/cond/switch_f

>layerfilter1_self_att_conv_head_0/bn/cond/control_dependency_1Identity2layerfilter1_self_att_conv_head_0/bn/cond/switch_f/^layerfilter1_self_att_conv_head_0/bn/cond/NoOp*
T0
*E
_class;
97loc:@layerfilter1_self_att_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
â
/layerfilter1_self_att_conv_head_0/bn/cond/MergeMerge>layerfilter1_self_att_conv_head_0/bn/cond/control_dependency_1<layerfilter1_self_att_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter1_self_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

4layerfilter1_self_att_conv_head_0/bn/cond_1/switch_tIdentity4layerfilter1_self_att_conv_head_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

4layerfilter1_self_att_conv_head_0/bn/cond_1/switch_fIdentity2layerfilter1_self_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
o
3layerfilter1_self_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ö
4layerfilter1_self_att_conv_head_0/bn/cond_1/IdentityIdentity=layerfilter1_self_att_conv_head_0/bn/cond_1/Identity/Switch:10^layerfilter1_self_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
¤
;layerfilter1_self_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch4layerfilter1_self_att_conv_head_0/bn/moments/Squeeze3layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter1_self_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
Ú
6layerfilter1_self_att_conv_head_0/bn/cond_1/Identity_1Identity?layerfilter1_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:10^layerfilter1_self_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
ª
=layerfilter1_self_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch6layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_13layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::

4layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_1Switchwlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
£
4layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_2Switchylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
1layerfilter1_self_att_conv_head_0/bn/cond_1/MergeMerge4layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_14layerfilter1_self_att_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
Ú
3layerfilter1_self_att_conv_head_0/bn/cond_1/Merge_1Merge4layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_26layerfilter1_self_att_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
y
4layerfilter1_self_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
É
2layerfilter1_self_att_conv_head_0/bn/batchnorm/addAdd3layerfilter1_self_att_conv_head_0/bn/cond_1/Merge_14layerfilter1_self_att_conv_head_0/bn/batchnorm/add/y*
T0*
_output_shapes
:

4layerfilter1_self_att_conv_head_0/bn/batchnorm/RsqrtRsqrt2layerfilter1_self_att_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
:
Å
2layerfilter1_self_att_conv_head_0/bn/batchnorm/mulMul4layerfilter1_self_att_conv_head_0/bn/batchnorm/Rsqrt/layerfilter1_self_att_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:
Ë
4layerfilter1_self_att_conv_head_0/bn/batchnorm/mul_1Mul)layerfilter1_self_att_conv_head_0/BiasAdd2layerfilter1_self_att_conv_head_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2
Ç
4layerfilter1_self_att_conv_head_0/bn/batchnorm/mul_2Mul1layerfilter1_self_att_conv_head_0/bn/cond_1/Merge2layerfilter1_self_att_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:
Ä
2layerfilter1_self_att_conv_head_0/bn/batchnorm/subSub.layerfilter1_self_att_conv_head_0/bn/beta/read4layerfilter1_self_att_conv_head_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:
Ö
4layerfilter1_self_att_conv_head_0/bn/batchnorm/add_1Add4layerfilter1_self_att_conv_head_0/bn/batchnorm/mul_12layerfilter1_self_att_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2

&layerfilter1_self_att_conv_head_0/ReluRelu4layerfilter1_self_att_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2
á
Jlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/shapeConst*%
valueB"      @      *
dtype0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*
_output_shapes
:
Ë
Hlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/minConst*
valueB
 *¾*
dtype0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*
_output_shapes
: 
Ë
Hlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/maxConst*
valueB
 *>*
dtype0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*
_output_shapes
: 
Ä
Rlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@
Â
Hlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/subSubHlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/maxHlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*
_output_shapes
: 
Ü
Hlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/mulMulRlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/RandomUniformHlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/sub*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@
Î
Dlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniformAddHlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/mulHlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@
ú
)layerfilter1_neib_att_conv_head_0/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@
Ò
0layerfilter1_neib_att_conv_head_0/weights/AssignAssign)layerfilter1_neib_att_conv_head_0/weightsDlayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@
ã
.layerfilter1_neib_att_conv_head_0/weights/readIdentity)layerfilter1_neib_att_conv_head_0/weights"/device:CPU:0*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@

(layerfilter1_neib_att_conv_head_0/L2LossL2Loss.layerfilter1_neib_att_conv_head_0/weights/read*
T0*
_output_shapes
: 
t
/layerfilter1_neib_att_conv_head_0/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter1_neib_att_conv_head_0/weight_lossMul(layerfilter1_neib_att_conv_head_0/L2Loss/layerfilter1_neib_att_conv_head_0/weight_loss/y*
T0*
_output_shapes
: 

(layerfilter1_neib_att_conv_head_0/Conv2DConv2Dlayerfilter1_edgefea_0/Relu.layerfilter1_neib_att_conv_head_0/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2

Ä
:layerfilter1_neib_att_conv_head_0/biases/Initializer/ConstConst*
valueB*    *
dtype0*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_0/biases*
_output_shapes
:
à
(layerfilter1_neib_att_conv_head_0/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *;
_class1
/-loc:@layerfilter1_neib_att_conv_head_0/biases*
_output_shapes
:
¹
/layerfilter1_neib_att_conv_head_0/biases/AssignAssign(layerfilter1_neib_att_conv_head_0/biases:layerfilter1_neib_att_conv_head_0/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_0/biases*
_output_shapes
:
Ô
-layerfilter1_neib_att_conv_head_0/biases/readIdentity(layerfilter1_neib_att_conv_head_0/biases"/device:CPU:0*
T0*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_0/biases*
_output_shapes
:
Õ
)layerfilter1_neib_att_conv_head_0/BiasAddBiasAdd(layerfilter1_neib_att_conv_head_0/Conv2D-layerfilter1_neib_att_conv_head_0/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2

w
*layerfilter1_neib_att_conv_head_0/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

)layerfilter1_neib_att_conv_head_0/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:

0layerfilter1_neib_att_conv_head_0/bn/beta/AssignAssign)layerfilter1_neib_att_conv_head_0/bn/beta*layerfilter1_neib_att_conv_head_0/bn/Const*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/bn/beta*
_output_shapes
:
È
.layerfilter1_neib_att_conv_head_0/bn/beta/readIdentity)layerfilter1_neib_att_conv_head_0/bn/beta*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/bn/beta*
_output_shapes
:
y
,layerfilter1_neib_att_conv_head_0/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

*layerfilter1_neib_att_conv_head_0/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
¢
1layerfilter1_neib_att_conv_head_0/bn/gamma/AssignAssign*layerfilter1_neib_att_conv_head_0/bn/gamma,layerfilter1_neib_att_conv_head_0/bn/Const_1*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_0/bn/gamma*
_output_shapes
:
Ë
/layerfilter1_neib_att_conv_head_0/bn/gamma/readIdentity*layerfilter1_neib_att_conv_head_0/bn/gamma*
T0*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_0/bn/gamma*
_output_shapes
:

Clayerfilter1_neib_att_conv_head_0/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
÷
1layerfilter1_neib_att_conv_head_0/bn/moments/meanMean)layerfilter1_neib_att_conv_head_0/BiasAddClayerfilter1_neib_att_conv_head_0/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
­
9layerfilter1_neib_att_conv_head_0/bn/moments/StopGradientStopGradient1layerfilter1_neib_att_conv_head_0/bn/moments/mean*
T0*&
_output_shapes
:
ê
>layerfilter1_neib_att_conv_head_0/bn/moments/SquaredDifferenceSquaredDifference)layerfilter1_neib_att_conv_head_0/BiasAdd9layerfilter1_neib_att_conv_head_0/bn/moments/StopGradient*
T0*&
_output_shapes
:2


Glayerfilter1_neib_att_conv_head_0/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

5layerfilter1_neib_att_conv_head_0/bn/moments/varianceMean>layerfilter1_neib_att_conv_head_0/bn/moments/SquaredDifferenceGlayerfilter1_neib_att_conv_head_0/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
°
4layerfilter1_neib_att_conv_head_0/bn/moments/SqueezeSqueeze1layerfilter1_neib_att_conv_head_0/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:
¶
6layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1Squeeze5layerfilter1_neib_att_conv_head_0/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:
{
0layerfilter1_neib_att_conv_head_0/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter1_neib_att_conv_head_0/bn/cond/switch_tIdentity2layerfilter1_neib_att_conv_head_0/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter1_neib_att_conv_head_0/bn/cond/switch_fIdentity0layerfilter1_neib_att_conv_head_0/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter1_neib_att_conv_head_0/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ç
layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ø
layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ð
layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
æ
rlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
ylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¤
wlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ë
layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ê
tlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
{layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ª
ylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Â
Hlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter1_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ú
Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter1_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ª
Vlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ô
Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
È
_layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter1_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
È
alayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze1layerfilter1_neib_att_conv_head_0/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
¼
Vlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
Rlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
À
Ylayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter1_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Þ
Zlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter1_neib_att_conv_head_0/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
°
Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ü
Zlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Î
alayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter1_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Î
clayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_11layerfilter1_neib_att_conv_head_0/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::
Ä
Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
Tlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
[layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter1_neib_att_conv_head_0/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
«
Blayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverageNoOpS^layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_13^layerfilter1_neib_att_conv_head_0/bn/cond/switch_t
©
<layerfilter1_neib_att_conv_head_0/bn/cond/control_dependencyIdentity2layerfilter1_neib_att_conv_head_0/bn/cond/switch_tC^layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter1_neib_att_conv_head_0/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter1_neib_att_conv_head_0/bn/cond/NoOpNoOp3^layerfilter1_neib_att_conv_head_0/bn/cond/switch_f

>layerfilter1_neib_att_conv_head_0/bn/cond/control_dependency_1Identity2layerfilter1_neib_att_conv_head_0/bn/cond/switch_f/^layerfilter1_neib_att_conv_head_0/bn/cond/NoOp*
T0
*E
_class;
97loc:@layerfilter1_neib_att_conv_head_0/bn/cond/switch_f*
_output_shapes
: 
â
/layerfilter1_neib_att_conv_head_0/bn/cond/MergeMerge>layerfilter1_neib_att_conv_head_0/bn/cond/control_dependency_1<layerfilter1_neib_att_conv_head_0/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter1_neib_att_conv_head_0/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

4layerfilter1_neib_att_conv_head_0/bn/cond_1/switch_tIdentity4layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

4layerfilter1_neib_att_conv_head_0/bn/cond_1/switch_fIdentity2layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch*
T0
*
_output_shapes
: 
o
3layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ö
4layerfilter1_neib_att_conv_head_0/bn/cond_1/IdentityIdentity=layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity/Switch:10^layerfilter1_neib_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
¤
;layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity/SwitchSwitch4layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze3layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze* 
_output_shapes
::
Ú
6layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity_1Identity?layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:10^layerfilter1_neib_att_conv_head_0/bn/cond/Merge*
T0*
_output_shapes
:
ª
=layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity_1/SwitchSwitch6layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_13layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1* 
_output_shapes
::

4layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_1Switchwlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
£
4layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_2Switchylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
1layerfilter1_neib_att_conv_head_0/bn/cond_1/MergeMerge4layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_14layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
Ú
3layerfilter1_neib_att_conv_head_0/bn/cond_1/Merge_1Merge4layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_26layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
y
4layerfilter1_neib_att_conv_head_0/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
É
2layerfilter1_neib_att_conv_head_0/bn/batchnorm/addAdd3layerfilter1_neib_att_conv_head_0/bn/cond_1/Merge_14layerfilter1_neib_att_conv_head_0/bn/batchnorm/add/y*
T0*
_output_shapes
:

4layerfilter1_neib_att_conv_head_0/bn/batchnorm/RsqrtRsqrt2layerfilter1_neib_att_conv_head_0/bn/batchnorm/add*
T0*
_output_shapes
:
Å
2layerfilter1_neib_att_conv_head_0/bn/batchnorm/mulMul4layerfilter1_neib_att_conv_head_0/bn/batchnorm/Rsqrt/layerfilter1_neib_att_conv_head_0/bn/gamma/read*
T0*
_output_shapes
:
Ë
4layerfilter1_neib_att_conv_head_0/bn/batchnorm/mul_1Mul)layerfilter1_neib_att_conv_head_0/BiasAdd2layerfilter1_neib_att_conv_head_0/bn/batchnorm/mul*
T0*&
_output_shapes
:2

Ç
4layerfilter1_neib_att_conv_head_0/bn/batchnorm/mul_2Mul1layerfilter1_neib_att_conv_head_0/bn/cond_1/Merge2layerfilter1_neib_att_conv_head_0/bn/batchnorm/mul*
T0*
_output_shapes
:
Ä
2layerfilter1_neib_att_conv_head_0/bn/batchnorm/subSub.layerfilter1_neib_att_conv_head_0/bn/beta/read4layerfilter1_neib_att_conv_head_0/bn/batchnorm/mul_2*
T0*
_output_shapes
:
Ö
4layerfilter1_neib_att_conv_head_0/bn/batchnorm/add_1Add4layerfilter1_neib_att_conv_head_0/bn/batchnorm/mul_12layerfilter1_neib_att_conv_head_0/bn/batchnorm/sub*
T0*&
_output_shapes
:2


&layerfilter1_neib_att_conv_head_0/ReluRelu4layerfilter1_neib_att_conv_head_0/bn/batchnorm/add_1*
T0*&
_output_shapes
:2


add_10Add&layerfilter1_self_att_conv_head_0/Relu&layerfilter1_neib_att_conv_head_0/Relu*
T0*&
_output_shapes
:2

i
transpose_7/permConst*%
valueB"             *
dtype0*
_output_shapes
:
p
transpose_7	Transposeadd_10transpose_7/perm*
T0*
Tperm0*&
_output_shapes
:2

V
LeakyRelu_1/alphaConst*
valueB
 *ÍÌL>*
dtype0*
_output_shapes
: 
g
LeakyRelu_1/mulMulLeakyRelu_1/alphatranspose_7*
T0*&
_output_shapes
:2

m
LeakyRelu_1/MaximumMaximumLeakyRelu_1/multranspose_7*
T0*&
_output_shapes
:2

`
Shape_3Const*%
valueB"   2      
   *
dtype0*
_output_shapes
:
H
Rank_1Const*
value	B :*
dtype0*
_output_shapes
: 
`
Shape_4Const*%
valueB"   2      
   *
dtype0*
_output_shapes
:
I
Sub_1/yConst*
value	B :*
dtype0*
_output_shapes
: 
>
Sub_1SubRank_1Sub_1/y*
T0*
_output_shapes
: 
V
Slice_1/beginPackSub_1*
N*
T0*

axis *
_output_shapes
:
V
Slice_1/sizeConst*
valueB:*
dtype0*
_output_shapes
:
h
Slice_1SliceShape_4Slice_1/beginSlice_1/size*
T0*
Index0*
_output_shapes
:
d
concat_4/values_0Const*
valueB:
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
:
O
concat_4/axisConst*
value	B : *
dtype0*
_output_shapes
: 
y
concat_4ConcatV2concat_4/values_0Slice_1concat_4/axis*
N*
T0*

Tidx0*
_output_shapes
:
j
	Reshape_6ReshapeLeakyRelu_1/Maximumconcat_4*
T0*
Tshape0*
_output_shapes

:2

H
	Softmax_1Softmax	Reshape_6*
T0*
_output_shapes

:2

g
	Reshape_7Reshape	Softmax_1Shape_3*
T0*
Tshape0*&
_output_shapes
:2


MatMul_3BatchMatMul	Reshape_7layerfilter1_edgefea_0/Relu*
T0*
adj_x( *
adj_y( *&
_output_shapes
:2@
¡
2BiasAdd_1/biases/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*#
_class
loc:@BiasAdd_1/biases*
_output_shapes
:

(BiasAdd_1/biases/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*#
_class
loc:@BiasAdd_1/biases*
_output_shapes
: 
ä
"BiasAdd_1/biases/Initializer/zerosFill2BiasAdd_1/biases/Initializer/zeros/shape_as_tensor(BiasAdd_1/biases/Initializer/zeros/Const*
T0*

index_type0*#
_class
loc:@BiasAdd_1/biases*
_output_shapes
:@
¡
BiasAdd_1/biases
VariableV2*
shape:@*
dtype0*
	container *
shared_name *#
_class
loc:@BiasAdd_1/biases*
_output_shapes
:@
Ê
BiasAdd_1/biases/AssignAssignBiasAdd_1/biases"BiasAdd_1/biases/Initializer/zeros*
T0*
validate_shape(*
use_locking(*#
_class
loc:@BiasAdd_1/biases*
_output_shapes
:@
}
BiasAdd_1/biases/readIdentityBiasAdd_1/biases*
T0*#
_class
loc:@BiasAdd_1/biases*
_output_shapes
:@

BiasAdd_1/BiasAddBiasAddMatMul_3BiasAdd_1/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2@
R
Relu_1ReluBiasAdd_1/BiasAdd*
T0*&
_output_shapes
:2@
a
	Squeeze_6Squeezegapnet01/Relu*
T0*
squeeze_dims
 *
_output_shapes
:	2
S
ExpandDims_12/dimConst*
value	B : *
dtype0*
_output_shapes
: 
s
ExpandDims_12
ExpandDims	Squeeze_6ExpandDims_12/dim*
T0*

Tdim0*#
_output_shapes
:2
\
ExpandDims_13/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
{
ExpandDims_13
ExpandDimsExpandDims_12ExpandDims_13/dim*
T0*

Tdim0*'
_output_shapes
:2
Ý
Hlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/shapeConst*%
valueB"         @   *
dtype0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*
_output_shapes
:
Ç
Flayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/minConst*
valueB
 *ó5¾*
dtype0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*
_output_shapes
: 
Ç
Flayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/maxConst*
valueB
 *ó5>*
dtype0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*
_output_shapes
: 
¿
Playerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/RandomUniformRandomUniformHlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@
º
Flayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/subSubFlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/maxFlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/min*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*
_output_shapes
: 
Õ
Flayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/mulMulPlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/RandomUniformFlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/sub*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@
Ç
Blayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniformAddFlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/mulFlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform/min*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@
ø
'layerfilter1_newfea_conv_head_1/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@
Ë
.layerfilter1_newfea_conv_head_1/weights/AssignAssign'layerfilter1_newfea_conv_head_1/weightsBlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@
Þ
,layerfilter1_newfea_conv_head_1/weights/readIdentity'layerfilter1_newfea_conv_head_1/weights"/device:CPU:0*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@

&layerfilter1_newfea_conv_head_1/L2LossL2Loss,layerfilter1_newfea_conv_head_1/weights/read*
T0*
_output_shapes
: 
r
-layerfilter1_newfea_conv_head_1/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
ª
+layerfilter1_newfea_conv_head_1/weight_lossMul&layerfilter1_newfea_conv_head_1/L2Loss-layerfilter1_newfea_conv_head_1/weight_loss/y*
T0*
_output_shapes
: 

&layerfilter1_newfea_conv_head_1/Conv2DConv2DExpandDims_13,layerfilter1_newfea_conv_head_1/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2@
u
(layerfilter1_newfea_conv_head_1/bn/ConstConst*
valueB@*    *
dtype0*
_output_shapes
:@

'layerfilter1_newfea_conv_head_1/bn/beta
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@

.layerfilter1_newfea_conv_head_1/bn/beta/AssignAssign'layerfilter1_newfea_conv_head_1/bn/beta(layerfilter1_newfea_conv_head_1/bn/Const*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/bn/beta*
_output_shapes
:@
Â
,layerfilter1_newfea_conv_head_1/bn/beta/readIdentity'layerfilter1_newfea_conv_head_1/bn/beta*
T0*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/bn/beta*
_output_shapes
:@
w
*layerfilter1_newfea_conv_head_1/bn/Const_1Const*
valueB@*  ?*
dtype0*
_output_shapes
:@

(layerfilter1_newfea_conv_head_1/bn/gamma
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@

/layerfilter1_newfea_conv_head_1/bn/gamma/AssignAssign(layerfilter1_newfea_conv_head_1/bn/gamma*layerfilter1_newfea_conv_head_1/bn/Const_1*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_newfea_conv_head_1/bn/gamma*
_output_shapes
:@
Å
-layerfilter1_newfea_conv_head_1/bn/gamma/readIdentity(layerfilter1_newfea_conv_head_1/bn/gamma*
T0*;
_class1
/-loc:@layerfilter1_newfea_conv_head_1/bn/gamma*
_output_shapes
:@

Alayerfilter1_newfea_conv_head_1/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ð
/layerfilter1_newfea_conv_head_1/bn/moments/meanMean&layerfilter1_newfea_conv_head_1/Conv2DAlayerfilter1_newfea_conv_head_1/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@
©
7layerfilter1_newfea_conv_head_1/bn/moments/StopGradientStopGradient/layerfilter1_newfea_conv_head_1/bn/moments/mean*
T0*&
_output_shapes
:@
ã
<layerfilter1_newfea_conv_head_1/bn/moments/SquaredDifferenceSquaredDifference&layerfilter1_newfea_conv_head_1/Conv2D7layerfilter1_newfea_conv_head_1/bn/moments/StopGradient*
T0*&
_output_shapes
:2@

Elayerfilter1_newfea_conv_head_1/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

3layerfilter1_newfea_conv_head_1/bn/moments/varianceMean<layerfilter1_newfea_conv_head_1/bn/moments/SquaredDifferenceElayerfilter1_newfea_conv_head_1/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@
¬
2layerfilter1_newfea_conv_head_1/bn/moments/SqueezeSqueeze/layerfilter1_newfea_conv_head_1/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:@
²
4layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1Squeeze3layerfilter1_newfea_conv_head_1/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:@
y
.layerfilter1_newfea_conv_head_1/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

0layerfilter1_newfea_conv_head_1/bn/cond/switch_tIdentity0layerfilter1_newfea_conv_head_1/bn/cond/Switch:1*
T0
*
_output_shapes
: 

0layerfilter1_newfea_conv_head_1/bn/cond/switch_fIdentity.layerfilter1_newfea_conv_head_1/bn/cond/Switch*
T0
*
_output_shapes
: 
k
/layerfilter1_newfea_conv_head_1/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ß
layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ð
layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
à
layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Þ
nlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ä
ulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

slayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/readIdentitynlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ã
layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ô
layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
è
layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
â
playerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Ì
wlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

ulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
¾
Flayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/decayConst1^layerfilter1_newfea_conv_head_1/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ò
Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst1^layerfilter1_newfea_conv_head_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
 
Tlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubVlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xFlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/decay*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ê
Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Sub_layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1alayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
¼
]layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchslayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read/layerfilter1_newfea_conv_head_1/bn/cond/pred_id*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
À
_layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch2layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/layerfilter1_newfea_conv_head_1/bn/cond/pred_id*
T0*E
_class;
97loc:@layerfilter1_newfea_conv_head_1/bn/moments/Squeeze* 
_output_shapes
:@:@
²
Tlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulVlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Tlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ê
Playerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubYlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Tlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
´
Wlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/layerfilter1_newfea_conv_head_1/bn/cond/pred_id*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
Ö
Xlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst1^layerfilter1_newfea_conv_head_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¦
Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubXlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xFlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/decay*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ò
Xlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subalayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1clayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Â
_layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read/layerfilter1_newfea_conv_head_1/bn/cond/pred_id*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Æ
alayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch4layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/layerfilter1_newfea_conv_head_1/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1* 
_output_shapes
:@:@
º
Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulXlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Ò
Rlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub[layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
º
Ylayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/layerfilter1_newfea_conv_head_1/bn/cond/pred_id*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
£
@layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverageNoOpQ^layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvgS^layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_11^layerfilter1_newfea_conv_head_1/bn/cond/switch_t
¡
:layerfilter1_newfea_conv_head_1/bn/cond/control_dependencyIdentity0layerfilter1_newfea_conv_head_1/bn/cond/switch_tA^layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage*
T0
*C
_class9
75loc:@layerfilter1_newfea_conv_head_1/bn/cond/switch_t*
_output_shapes
: 
g
,layerfilter1_newfea_conv_head_1/bn/cond/NoOpNoOp1^layerfilter1_newfea_conv_head_1/bn/cond/switch_f

<layerfilter1_newfea_conv_head_1/bn/cond/control_dependency_1Identity0layerfilter1_newfea_conv_head_1/bn/cond/switch_f-^layerfilter1_newfea_conv_head_1/bn/cond/NoOp*
T0
*C
_class9
75loc:@layerfilter1_newfea_conv_head_1/bn/cond/switch_f*
_output_shapes
: 
Ü
-layerfilter1_newfea_conv_head_1/bn/cond/MergeMerge<layerfilter1_newfea_conv_head_1/bn/cond/control_dependency_1:layerfilter1_newfea_conv_head_1/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
{
0layerfilter1_newfea_conv_head_1/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter1_newfea_conv_head_1/bn/cond_1/switch_tIdentity2layerfilter1_newfea_conv_head_1/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

2layerfilter1_newfea_conv_head_1/bn/cond_1/switch_fIdentity0layerfilter1_newfea_conv_head_1/bn/cond_1/Switch*
T0
*
_output_shapes
: 
m
1layerfilter1_newfea_conv_head_1/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ð
2layerfilter1_newfea_conv_head_1/bn/cond_1/IdentityIdentity;layerfilter1_newfea_conv_head_1/bn/cond_1/Identity/Switch:1.^layerfilter1_newfea_conv_head_1/bn/cond/Merge*
T0*
_output_shapes
:@

9layerfilter1_newfea_conv_head_1/bn/cond_1/Identity/SwitchSwitch2layerfilter1_newfea_conv_head_1/bn/moments/Squeeze1layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id*
T0*E
_class;
97loc:@layerfilter1_newfea_conv_head_1/bn/moments/Squeeze* 
_output_shapes
:@:@
Ô
4layerfilter1_newfea_conv_head_1/bn/cond_1/Identity_1Identity=layerfilter1_newfea_conv_head_1/bn/cond_1/Identity_1/Switch:1.^layerfilter1_newfea_conv_head_1/bn/cond/Merge*
T0*
_output_shapes
:@
¢
;layerfilter1_newfea_conv_head_1/bn/cond_1/Identity_1/SwitchSwitch4layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_11layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1* 
_output_shapes
:@:@

2layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_1Switchslayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id*
T0*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@

2layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_2Switchulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id*
T0*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
Ð
/layerfilter1_newfea_conv_head_1/bn/cond_1/MergeMerge2layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_12layerfilter1_newfea_conv_head_1/bn/cond_1/Identity*
T0*
N*
_output_shapes

:@: 
Ô
1layerfilter1_newfea_conv_head_1/bn/cond_1/Merge_1Merge2layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_24layerfilter1_newfea_conv_head_1/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:@: 
w
2layerfilter1_newfea_conv_head_1/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
Ã
0layerfilter1_newfea_conv_head_1/bn/batchnorm/addAdd1layerfilter1_newfea_conv_head_1/bn/cond_1/Merge_12layerfilter1_newfea_conv_head_1/bn/batchnorm/add/y*
T0*
_output_shapes
:@

2layerfilter1_newfea_conv_head_1/bn/batchnorm/RsqrtRsqrt0layerfilter1_newfea_conv_head_1/bn/batchnorm/add*
T0*
_output_shapes
:@
¿
0layerfilter1_newfea_conv_head_1/bn/batchnorm/mulMul2layerfilter1_newfea_conv_head_1/bn/batchnorm/Rsqrt-layerfilter1_newfea_conv_head_1/bn/gamma/read*
T0*
_output_shapes
:@
Ä
2layerfilter1_newfea_conv_head_1/bn/batchnorm/mul_1Mul&layerfilter1_newfea_conv_head_1/Conv2D0layerfilter1_newfea_conv_head_1/bn/batchnorm/mul*
T0*&
_output_shapes
:2@
Á
2layerfilter1_newfea_conv_head_1/bn/batchnorm/mul_2Mul/layerfilter1_newfea_conv_head_1/bn/cond_1/Merge0layerfilter1_newfea_conv_head_1/bn/batchnorm/mul*
T0*
_output_shapes
:@
¾
0layerfilter1_newfea_conv_head_1/bn/batchnorm/subSub,layerfilter1_newfea_conv_head_1/bn/beta/read2layerfilter1_newfea_conv_head_1/bn/batchnorm/mul_2*
T0*
_output_shapes
:@
Ð
2layerfilter1_newfea_conv_head_1/bn/batchnorm/add_1Add2layerfilter1_newfea_conv_head_1/bn/batchnorm/mul_10layerfilter1_newfea_conv_head_1/bn/batchnorm/sub*
T0*&
_output_shapes
:2@

$layerfilter1_newfea_conv_head_1/ReluRelu2layerfilter1_newfea_conv_head_1/bn/batchnorm/add_1*
T0*&
_output_shapes
:2@
a
	Squeeze_7SqueezeExpandDims_13*
T0*
squeeze_dims
 *
_output_shapes
:	2
S
ExpandDims_14/dimConst*
value	B : *
dtype0*
_output_shapes
: 
s
ExpandDims_14
ExpandDims	Squeeze_7ExpandDims_14/dim*
T0*

Tdim0*#
_output_shapes
:2
O
range_2/startConst*
value	B : *
dtype0*
_output_shapes
: 
O
range_2/limitConst*
value	B :*
dtype0*
_output_shapes
: 
O
range_2/deltaConst*
value	B :*
dtype0*
_output_shapes
: 
e
range_2Rangerange_2/startrange_2/limitrange_2/delta*

Tidx0*
_output_shapes
:
I
mul_6/yConst*
value	B :2*
dtype0*
_output_shapes
: 
C
mul_6Mulrange_2mul_6/y*
T0*
_output_shapes
:
d
Reshape_8/shapeConst*!
valueB"         *
dtype0*
_output_shapes
:
g
	Reshape_8Reshapemul_6Reshape_8/shape*
T0*
Tshape0*"
_output_shapes
:
`
Reshape_9/shapeConst*
valueB"ÿÿÿÿ   *
dtype0*
_output_shapes
:
l
	Reshape_9ReshapeExpandDims_14Reshape_9/shape*
T0*
Tshape0*
_output_shapes
:	2
Q
add_11Add
TopKV2_1:1	Reshape_8*
T0*"
_output_shapes
:2


Gather_2Gather	Reshape_9add_11*
validate_indices(*
Tparams0*
Tindices0*'
_output_shapes
:2

i
Tile_3/multiplesConst*%
valueB"      
      *
dtype0*
_output_shapes
:
s
Tile_3TileExpandDims_13Tile_3/multiples*
T0*

Tmultiples0*'
_output_shapes
:2

P
sub_5SubTile_3Gather_2*
T0*'
_output_shapes
:2

Ë
?layerfilter1_edgefea_1/weights/Initializer/random_uniform/shapeConst*%
valueB"         @   *
dtype0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*
_output_shapes
:
µ
=layerfilter1_edgefea_1/weights/Initializer/random_uniform/minConst*
valueB
 *ó5¾*
dtype0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*
_output_shapes
: 
µ
=layerfilter1_edgefea_1/weights/Initializer/random_uniform/maxConst*
valueB
 *ó5>*
dtype0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*
_output_shapes
: 
¤
Glayerfilter1_edgefea_1/weights/Initializer/random_uniform/RandomUniformRandomUniform?layerfilter1_edgefea_1/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@

=layerfilter1_edgefea_1/weights/Initializer/random_uniform/subSub=layerfilter1_edgefea_1/weights/Initializer/random_uniform/max=layerfilter1_edgefea_1/weights/Initializer/random_uniform/min*
T0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*
_output_shapes
: 
±
=layerfilter1_edgefea_1/weights/Initializer/random_uniform/mulMulGlayerfilter1_edgefea_1/weights/Initializer/random_uniform/RandomUniform=layerfilter1_edgefea_1/weights/Initializer/random_uniform/sub*
T0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@
£
9layerfilter1_edgefea_1/weights/Initializer/random_uniformAdd=layerfilter1_edgefea_1/weights/Initializer/random_uniform/mul=layerfilter1_edgefea_1/weights/Initializer/random_uniform/min*
T0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@
æ
layerfilter1_edgefea_1/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@
§
%layerfilter1_edgefea_1/weights/AssignAssignlayerfilter1_edgefea_1/weights9layerfilter1_edgefea_1/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@
Ã
#layerfilter1_edgefea_1/weights/readIdentitylayerfilter1_edgefea_1/weights"/device:CPU:0*
T0*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@
m
layerfilter1_edgefea_1/L2LossL2Loss#layerfilter1_edgefea_1/weights/read*
T0*
_output_shapes
: 
i
$layerfilter1_edgefea_1/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 

"layerfilter1_edgefea_1/weight_lossMullayerfilter1_edgefea_1/L2Loss$layerfilter1_edgefea_1/weight_loss/y*
T0*
_output_shapes
: 
ó
layerfilter1_edgefea_1/Conv2DConv2Dsub_5#layerfilter1_edgefea_1/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2
@
®
/layerfilter1_edgefea_1/biases/Initializer/ConstConst*
valueB@*    *
dtype0*0
_class&
$"loc:@layerfilter1_edgefea_1/biases*
_output_shapes
:@
Ê
layerfilter1_edgefea_1/biases
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *0
_class&
$"loc:@layerfilter1_edgefea_1/biases*
_output_shapes
:@

$layerfilter1_edgefea_1/biases/AssignAssignlayerfilter1_edgefea_1/biases/layerfilter1_edgefea_1/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter1_edgefea_1/biases*
_output_shapes
:@
³
"layerfilter1_edgefea_1/biases/readIdentitylayerfilter1_edgefea_1/biases"/device:CPU:0*
T0*0
_class&
$"loc:@layerfilter1_edgefea_1/biases*
_output_shapes
:@
´
layerfilter1_edgefea_1/BiasAddBiasAddlayerfilter1_edgefea_1/Conv2D"layerfilter1_edgefea_1/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
@
l
layerfilter1_edgefea_1/bn/ConstConst*
valueB@*    *
dtype0*
_output_shapes
:@

layerfilter1_edgefea_1/bn/beta
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@
ñ
%layerfilter1_edgefea_1/bn/beta/AssignAssignlayerfilter1_edgefea_1/bn/betalayerfilter1_edgefea_1/bn/Const*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_1/bn/beta*
_output_shapes
:@
§
#layerfilter1_edgefea_1/bn/beta/readIdentitylayerfilter1_edgefea_1/bn/beta*
T0*1
_class'
%#loc:@layerfilter1_edgefea_1/bn/beta*
_output_shapes
:@
n
!layerfilter1_edgefea_1/bn/Const_1Const*
valueB@*  ?*
dtype0*
_output_shapes
:@

layerfilter1_edgefea_1/bn/gamma
VariableV2*
shape:@*
dtype0*
	container *
shared_name *
_output_shapes
:@
ö
&layerfilter1_edgefea_1/bn/gamma/AssignAssignlayerfilter1_edgefea_1/bn/gamma!layerfilter1_edgefea_1/bn/Const_1*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter1_edgefea_1/bn/gamma*
_output_shapes
:@
ª
$layerfilter1_edgefea_1/bn/gamma/readIdentitylayerfilter1_edgefea_1/bn/gamma*
T0*2
_class(
&$loc:@layerfilter1_edgefea_1/bn/gamma*
_output_shapes
:@

8layerfilter1_edgefea_1/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ö
&layerfilter1_edgefea_1/bn/moments/meanMeanlayerfilter1_edgefea_1/BiasAdd8layerfilter1_edgefea_1/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@

.layerfilter1_edgefea_1/bn/moments/StopGradientStopGradient&layerfilter1_edgefea_1/bn/moments/mean*
T0*&
_output_shapes
:@
É
3layerfilter1_edgefea_1/bn/moments/SquaredDifferenceSquaredDifferencelayerfilter1_edgefea_1/BiasAdd.layerfilter1_edgefea_1/bn/moments/StopGradient*
T0*&
_output_shapes
:2
@

<layerfilter1_edgefea_1/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
ó
*layerfilter1_edgefea_1/bn/moments/varianceMean3layerfilter1_edgefea_1/bn/moments/SquaredDifference<layerfilter1_edgefea_1/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:@

)layerfilter1_edgefea_1/bn/moments/SqueezeSqueeze&layerfilter1_edgefea_1/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:@
 
+layerfilter1_edgefea_1/bn/moments/Squeeze_1Squeeze*layerfilter1_edgefea_1/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:@
p
%layerfilter1_edgefea_1/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
}
'layerfilter1_edgefea_1/bn/cond/switch_tIdentity'layerfilter1_edgefea_1/bn/cond/Switch:1*
T0
*
_output_shapes
: 
{
'layerfilter1_edgefea_1/bn/cond/switch_fIdentity%layerfilter1_edgefea_1/bn/cond/Switch*
T0
*
_output_shapes
: 
b
&layerfilter1_edgefea_1/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
¹
~layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ª
tlayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

nlayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFill~layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensortlayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
¹
\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ú
clayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAveragenlayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
á
alayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
T0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
¾
layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
®
vlayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

playerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorvlayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
½
^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:@*
dtype0*
	container *
shared_name *q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

elayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssign^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverageplayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
ç
clayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentity^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
¬
=layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/decayConst(^layerfilter1_edgefea_1/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
­
Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst(^layerfilter1_edgefea_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ò
Klayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubMlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x=layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/decay*
T0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubVlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Xlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

Tlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchalayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/read&layerfilter1_edgefea_1/bn/cond/pred_id*
T0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@

Vlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch)layerfilter1_edgefea_1/bn/moments/Squeeze&layerfilter1_edgefea_1/bn/cond/pred_id*
T0*<
_class2
0.loc:@layerfilter1_edgefea_1/bn/moments/Squeeze* 
_output_shapes
:@:@

Klayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulMlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Klayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

Glayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubPlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Klayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ý
Nlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage&layerfilter1_edgefea_1/bn/cond/pred_id*
T0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
±
Olayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst(^layerfilter1_edgefea_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubOlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x=layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/decay*
T0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
¤
Olayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubXlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Zlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Vlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchclayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read&layerfilter1_edgefea_1/bn/cond/pred_id*
T0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
¢
Xlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch+layerfilter1_edgefea_1/bn/moments/Squeeze_1&layerfilter1_edgefea_1/bn/cond/pred_id*
T0*>
_class4
20loc:@layerfilter1_edgefea_1/bn/moments/Squeeze_1* 
_output_shapes
:@:@

Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulOlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
¤
Ilayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubRlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@

Playerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitch^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage&layerfilter1_edgefea_1/bn/cond/pred_id*
T0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
ÿ
7layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverageNoOpH^layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvgJ^layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1(^layerfilter1_edgefea_1/bn/cond/switch_t
ý
1layerfilter1_edgefea_1/bn/cond/control_dependencyIdentity'layerfilter1_edgefea_1/bn/cond/switch_t8^layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage*
T0
*:
_class0
.,loc:@layerfilter1_edgefea_1/bn/cond/switch_t*
_output_shapes
: 
U
#layerfilter1_edgefea_1/bn/cond/NoOpNoOp(^layerfilter1_edgefea_1/bn/cond/switch_f
ë
3layerfilter1_edgefea_1/bn/cond/control_dependency_1Identity'layerfilter1_edgefea_1/bn/cond/switch_f$^layerfilter1_edgefea_1/bn/cond/NoOp*
T0
*:
_class0
.,loc:@layerfilter1_edgefea_1/bn/cond/switch_f*
_output_shapes
: 
Á
$layerfilter1_edgefea_1/bn/cond/MergeMerge3layerfilter1_edgefea_1/bn/cond/control_dependency_11layerfilter1_edgefea_1/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
r
'layerfilter1_edgefea_1/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

)layerfilter1_edgefea_1/bn/cond_1/switch_tIdentity)layerfilter1_edgefea_1/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

)layerfilter1_edgefea_1/bn/cond_1/switch_fIdentity'layerfilter1_edgefea_1/bn/cond_1/Switch*
T0
*
_output_shapes
: 
d
(layerfilter1_edgefea_1/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
µ
)layerfilter1_edgefea_1/bn/cond_1/IdentityIdentity2layerfilter1_edgefea_1/bn/cond_1/Identity/Switch:1%^layerfilter1_edgefea_1/bn/cond/Merge*
T0*
_output_shapes
:@
ø
0layerfilter1_edgefea_1/bn/cond_1/Identity/SwitchSwitch)layerfilter1_edgefea_1/bn/moments/Squeeze(layerfilter1_edgefea_1/bn/cond_1/pred_id*
T0*<
_class2
0.loc:@layerfilter1_edgefea_1/bn/moments/Squeeze* 
_output_shapes
:@:@
¹
+layerfilter1_edgefea_1/bn/cond_1/Identity_1Identity4layerfilter1_edgefea_1/bn/cond_1/Identity_1/Switch:1%^layerfilter1_edgefea_1/bn/cond/Merge*
T0*
_output_shapes
:@
þ
2layerfilter1_edgefea_1/bn/cond_1/Identity_1/SwitchSwitch+layerfilter1_edgefea_1/bn/moments/Squeeze_1(layerfilter1_edgefea_1/bn/cond_1/pred_id*
T0*>
_class4
20loc:@layerfilter1_edgefea_1/bn/moments/Squeeze_1* 
_output_shapes
:@:@
Ü
)layerfilter1_edgefea_1/bn/cond_1/Switch_1Switchalayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/read(layerfilter1_edgefea_1/bn/cond_1/pred_id*
T0*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
:@:@
à
)layerfilter1_edgefea_1/bn/cond_1/Switch_2Switchclayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read(layerfilter1_edgefea_1/bn/cond_1/pred_id*
T0*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
:@:@
µ
&layerfilter1_edgefea_1/bn/cond_1/MergeMerge)layerfilter1_edgefea_1/bn/cond_1/Switch_1)layerfilter1_edgefea_1/bn/cond_1/Identity*
T0*
N*
_output_shapes

:@: 
¹
(layerfilter1_edgefea_1/bn/cond_1/Merge_1Merge)layerfilter1_edgefea_1/bn/cond_1/Switch_2+layerfilter1_edgefea_1/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:@: 
n
)layerfilter1_edgefea_1/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
¨
'layerfilter1_edgefea_1/bn/batchnorm/addAdd(layerfilter1_edgefea_1/bn/cond_1/Merge_1)layerfilter1_edgefea_1/bn/batchnorm/add/y*
T0*
_output_shapes
:@

)layerfilter1_edgefea_1/bn/batchnorm/RsqrtRsqrt'layerfilter1_edgefea_1/bn/batchnorm/add*
T0*
_output_shapes
:@
¤
'layerfilter1_edgefea_1/bn/batchnorm/mulMul)layerfilter1_edgefea_1/bn/batchnorm/Rsqrt$layerfilter1_edgefea_1/bn/gamma/read*
T0*
_output_shapes
:@
ª
)layerfilter1_edgefea_1/bn/batchnorm/mul_1Mullayerfilter1_edgefea_1/BiasAdd'layerfilter1_edgefea_1/bn/batchnorm/mul*
T0*&
_output_shapes
:2
@
¦
)layerfilter1_edgefea_1/bn/batchnorm/mul_2Mul&layerfilter1_edgefea_1/bn/cond_1/Merge'layerfilter1_edgefea_1/bn/batchnorm/mul*
T0*
_output_shapes
:@
£
'layerfilter1_edgefea_1/bn/batchnorm/subSub#layerfilter1_edgefea_1/bn/beta/read)layerfilter1_edgefea_1/bn/batchnorm/mul_2*
T0*
_output_shapes
:@
µ
)layerfilter1_edgefea_1/bn/batchnorm/add_1Add)layerfilter1_edgefea_1/bn/batchnorm/mul_1'layerfilter1_edgefea_1/bn/batchnorm/sub*
T0*&
_output_shapes
:2
@

layerfilter1_edgefea_1/ReluRelu)layerfilter1_edgefea_1/bn/batchnorm/add_1*
T0*&
_output_shapes
:2
@
á
Jlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/shapeConst*%
valueB"      @      *
dtype0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*
_output_shapes
:
Ë
Hlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/minConst*
valueB
 *¾*
dtype0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*
_output_shapes
: 
Ë
Hlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/maxConst*
valueB
 *>*
dtype0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*
_output_shapes
: 
Ä
Rlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@
Â
Hlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/subSubHlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/maxHlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*
_output_shapes
: 
Ü
Hlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/mulMulRlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/RandomUniformHlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/sub*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@
Î
Dlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniformAddHlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/mulHlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@
ú
)layerfilter1_self_att_conv_head_1/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@
Ò
0layerfilter1_self_att_conv_head_1/weights/AssignAssign)layerfilter1_self_att_conv_head_1/weightsDlayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@
ã
.layerfilter1_self_att_conv_head_1/weights/readIdentity)layerfilter1_self_att_conv_head_1/weights"/device:CPU:0*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@

(layerfilter1_self_att_conv_head_1/L2LossL2Loss.layerfilter1_self_att_conv_head_1/weights/read*
T0*
_output_shapes
: 
t
/layerfilter1_self_att_conv_head_1/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter1_self_att_conv_head_1/weight_lossMul(layerfilter1_self_att_conv_head_1/L2Loss/layerfilter1_self_att_conv_head_1/weight_loss/y*
T0*
_output_shapes
: 
¨
(layerfilter1_self_att_conv_head_1/Conv2DConv2D$layerfilter1_newfea_conv_head_1/Relu.layerfilter1_self_att_conv_head_1/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2
Ä
:layerfilter1_self_att_conv_head_1/biases/Initializer/ConstConst*
valueB*    *
dtype0*;
_class1
/-loc:@layerfilter1_self_att_conv_head_1/biases*
_output_shapes
:
à
(layerfilter1_self_att_conv_head_1/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *;
_class1
/-loc:@layerfilter1_self_att_conv_head_1/biases*
_output_shapes
:
¹
/layerfilter1_self_att_conv_head_1/biases/AssignAssign(layerfilter1_self_att_conv_head_1/biases:layerfilter1_self_att_conv_head_1/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_self_att_conv_head_1/biases*
_output_shapes
:
Ô
-layerfilter1_self_att_conv_head_1/biases/readIdentity(layerfilter1_self_att_conv_head_1/biases"/device:CPU:0*
T0*;
_class1
/-loc:@layerfilter1_self_att_conv_head_1/biases*
_output_shapes
:
Õ
)layerfilter1_self_att_conv_head_1/BiasAddBiasAdd(layerfilter1_self_att_conv_head_1/Conv2D-layerfilter1_self_att_conv_head_1/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
w
*layerfilter1_self_att_conv_head_1/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

)layerfilter1_self_att_conv_head_1/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:

0layerfilter1_self_att_conv_head_1/bn/beta/AssignAssign)layerfilter1_self_att_conv_head_1/bn/beta*layerfilter1_self_att_conv_head_1/bn/Const*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/bn/beta*
_output_shapes
:
È
.layerfilter1_self_att_conv_head_1/bn/beta/readIdentity)layerfilter1_self_att_conv_head_1/bn/beta*
T0*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/bn/beta*
_output_shapes
:
y
,layerfilter1_self_att_conv_head_1/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

*layerfilter1_self_att_conv_head_1/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
¢
1layerfilter1_self_att_conv_head_1/bn/gamma/AssignAssign*layerfilter1_self_att_conv_head_1/bn/gamma,layerfilter1_self_att_conv_head_1/bn/Const_1*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_self_att_conv_head_1/bn/gamma*
_output_shapes
:
Ë
/layerfilter1_self_att_conv_head_1/bn/gamma/readIdentity*layerfilter1_self_att_conv_head_1/bn/gamma*
T0*=
_class3
1/loc:@layerfilter1_self_att_conv_head_1/bn/gamma*
_output_shapes
:

Clayerfilter1_self_att_conv_head_1/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
÷
1layerfilter1_self_att_conv_head_1/bn/moments/meanMean)layerfilter1_self_att_conv_head_1/BiasAddClayerfilter1_self_att_conv_head_1/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
­
9layerfilter1_self_att_conv_head_1/bn/moments/StopGradientStopGradient1layerfilter1_self_att_conv_head_1/bn/moments/mean*
T0*&
_output_shapes
:
ê
>layerfilter1_self_att_conv_head_1/bn/moments/SquaredDifferenceSquaredDifference)layerfilter1_self_att_conv_head_1/BiasAdd9layerfilter1_self_att_conv_head_1/bn/moments/StopGradient*
T0*&
_output_shapes
:2

Glayerfilter1_self_att_conv_head_1/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

5layerfilter1_self_att_conv_head_1/bn/moments/varianceMean>layerfilter1_self_att_conv_head_1/bn/moments/SquaredDifferenceGlayerfilter1_self_att_conv_head_1/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
°
4layerfilter1_self_att_conv_head_1/bn/moments/SqueezeSqueeze1layerfilter1_self_att_conv_head_1/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:
¶
6layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1Squeeze5layerfilter1_self_att_conv_head_1/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:
{
0layerfilter1_self_att_conv_head_1/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter1_self_att_conv_head_1/bn/cond/switch_tIdentity2layerfilter1_self_att_conv_head_1/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter1_self_att_conv_head_1/bn/cond/switch_fIdentity0layerfilter1_self_att_conv_head_1/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter1_self_att_conv_head_1/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ç
layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ø
layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ð
layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
æ
rlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
ylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¤
wlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ë
layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ê
tlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
{layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ª
ylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Â
Hlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter1_self_att_conv_head_1/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ú
Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter1_self_att_conv_head_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ª
Vlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/decay*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ô
Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
È
_layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter1_self_att_conv_head_1/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
È
alayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter1_self_att_conv_head_1/bn/moments/Squeeze1layerfilter1_self_att_conv_head_1/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter1_self_att_conv_head_1/bn/moments/Squeeze* 
_output_shapes
::
¼
Vlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
Rlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
À
Ylayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter1_self_att_conv_head_1/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Þ
Zlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter1_self_att_conv_head_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
°
Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/decay*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ü
Zlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Î
alayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter1_self_att_conv_head_1/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Î
clayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_11layerfilter1_self_att_conv_head_1/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1* 
_output_shapes
::
Ä
Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
Tlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
[layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter1_self_att_conv_head_1/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
«
Blayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverageNoOpS^layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_13^layerfilter1_self_att_conv_head_1/bn/cond/switch_t
©
<layerfilter1_self_att_conv_head_1/bn/cond/control_dependencyIdentity2layerfilter1_self_att_conv_head_1/bn/cond/switch_tC^layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter1_self_att_conv_head_1/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter1_self_att_conv_head_1/bn/cond/NoOpNoOp3^layerfilter1_self_att_conv_head_1/bn/cond/switch_f

>layerfilter1_self_att_conv_head_1/bn/cond/control_dependency_1Identity2layerfilter1_self_att_conv_head_1/bn/cond/switch_f/^layerfilter1_self_att_conv_head_1/bn/cond/NoOp*
T0
*E
_class;
97loc:@layerfilter1_self_att_conv_head_1/bn/cond/switch_f*
_output_shapes
: 
â
/layerfilter1_self_att_conv_head_1/bn/cond/MergeMerge>layerfilter1_self_att_conv_head_1/bn/cond/control_dependency_1<layerfilter1_self_att_conv_head_1/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter1_self_att_conv_head_1/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

4layerfilter1_self_att_conv_head_1/bn/cond_1/switch_tIdentity4layerfilter1_self_att_conv_head_1/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

4layerfilter1_self_att_conv_head_1/bn/cond_1/switch_fIdentity2layerfilter1_self_att_conv_head_1/bn/cond_1/Switch*
T0
*
_output_shapes
: 
o
3layerfilter1_self_att_conv_head_1/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ö
4layerfilter1_self_att_conv_head_1/bn/cond_1/IdentityIdentity=layerfilter1_self_att_conv_head_1/bn/cond_1/Identity/Switch:10^layerfilter1_self_att_conv_head_1/bn/cond/Merge*
T0*
_output_shapes
:
¤
;layerfilter1_self_att_conv_head_1/bn/cond_1/Identity/SwitchSwitch4layerfilter1_self_att_conv_head_1/bn/moments/Squeeze3layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter1_self_att_conv_head_1/bn/moments/Squeeze* 
_output_shapes
::
Ú
6layerfilter1_self_att_conv_head_1/bn/cond_1/Identity_1Identity?layerfilter1_self_att_conv_head_1/bn/cond_1/Identity_1/Switch:10^layerfilter1_self_att_conv_head_1/bn/cond/Merge*
T0*
_output_shapes
:
ª
=layerfilter1_self_att_conv_head_1/bn/cond_1/Identity_1/SwitchSwitch6layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_13layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1* 
_output_shapes
::

4layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_1Switchwlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id*
T0*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
£
4layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_2Switchylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id*
T0*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
1layerfilter1_self_att_conv_head_1/bn/cond_1/MergeMerge4layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_14layerfilter1_self_att_conv_head_1/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
Ú
3layerfilter1_self_att_conv_head_1/bn/cond_1/Merge_1Merge4layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_26layerfilter1_self_att_conv_head_1/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
y
4layerfilter1_self_att_conv_head_1/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
É
2layerfilter1_self_att_conv_head_1/bn/batchnorm/addAdd3layerfilter1_self_att_conv_head_1/bn/cond_1/Merge_14layerfilter1_self_att_conv_head_1/bn/batchnorm/add/y*
T0*
_output_shapes
:

4layerfilter1_self_att_conv_head_1/bn/batchnorm/RsqrtRsqrt2layerfilter1_self_att_conv_head_1/bn/batchnorm/add*
T0*
_output_shapes
:
Å
2layerfilter1_self_att_conv_head_1/bn/batchnorm/mulMul4layerfilter1_self_att_conv_head_1/bn/batchnorm/Rsqrt/layerfilter1_self_att_conv_head_1/bn/gamma/read*
T0*
_output_shapes
:
Ë
4layerfilter1_self_att_conv_head_1/bn/batchnorm/mul_1Mul)layerfilter1_self_att_conv_head_1/BiasAdd2layerfilter1_self_att_conv_head_1/bn/batchnorm/mul*
T0*&
_output_shapes
:2
Ç
4layerfilter1_self_att_conv_head_1/bn/batchnorm/mul_2Mul1layerfilter1_self_att_conv_head_1/bn/cond_1/Merge2layerfilter1_self_att_conv_head_1/bn/batchnorm/mul*
T0*
_output_shapes
:
Ä
2layerfilter1_self_att_conv_head_1/bn/batchnorm/subSub.layerfilter1_self_att_conv_head_1/bn/beta/read4layerfilter1_self_att_conv_head_1/bn/batchnorm/mul_2*
T0*
_output_shapes
:
Ö
4layerfilter1_self_att_conv_head_1/bn/batchnorm/add_1Add4layerfilter1_self_att_conv_head_1/bn/batchnorm/mul_12layerfilter1_self_att_conv_head_1/bn/batchnorm/sub*
T0*&
_output_shapes
:2

&layerfilter1_self_att_conv_head_1/ReluRelu4layerfilter1_self_att_conv_head_1/bn/batchnorm/add_1*
T0*&
_output_shapes
:2
á
Jlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/shapeConst*%
valueB"      @      *
dtype0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*
_output_shapes
:
Ë
Hlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/minConst*
valueB
 *¾*
dtype0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*
_output_shapes
: 
Ë
Hlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/maxConst*
valueB
 *>*
dtype0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*
_output_shapes
: 
Ä
Rlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/RandomUniformRandomUniformJlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@
Â
Hlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/subSubHlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/maxHlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*
_output_shapes
: 
Ü
Hlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/mulMulRlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/RandomUniformHlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/sub*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@
Î
Dlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniformAddHlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/mulHlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform/min*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@
ú
)layerfilter1_neib_att_conv_head_1/weights
VariableV2"/device:CPU:0*
shape:@*
dtype0*
	container *
shared_name *<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@
Ò
0layerfilter1_neib_att_conv_head_1/weights/AssignAssign)layerfilter1_neib_att_conv_head_1/weightsDlayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@
ã
.layerfilter1_neib_att_conv_head_1/weights/readIdentity)layerfilter1_neib_att_conv_head_1/weights"/device:CPU:0*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@

(layerfilter1_neib_att_conv_head_1/L2LossL2Loss.layerfilter1_neib_att_conv_head_1/weights/read*
T0*
_output_shapes
: 
t
/layerfilter1_neib_att_conv_head_1/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
°
-layerfilter1_neib_att_conv_head_1/weight_lossMul(layerfilter1_neib_att_conv_head_1/L2Loss/layerfilter1_neib_att_conv_head_1/weight_loss/y*
T0*
_output_shapes
: 

(layerfilter1_neib_att_conv_head_1/Conv2DConv2Dlayerfilter1_edgefea_1/Relu.layerfilter1_neib_att_conv_head_1/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2

Ä
:layerfilter1_neib_att_conv_head_1/biases/Initializer/ConstConst*
valueB*    *
dtype0*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_1/biases*
_output_shapes
:
à
(layerfilter1_neib_att_conv_head_1/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *;
_class1
/-loc:@layerfilter1_neib_att_conv_head_1/biases*
_output_shapes
:
¹
/layerfilter1_neib_att_conv_head_1/biases/AssignAssign(layerfilter1_neib_att_conv_head_1/biases:layerfilter1_neib_att_conv_head_1/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_1/biases*
_output_shapes
:
Ô
-layerfilter1_neib_att_conv_head_1/biases/readIdentity(layerfilter1_neib_att_conv_head_1/biases"/device:CPU:0*
T0*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_1/biases*
_output_shapes
:
Õ
)layerfilter1_neib_att_conv_head_1/BiasAddBiasAdd(layerfilter1_neib_att_conv_head_1/Conv2D-layerfilter1_neib_att_conv_head_1/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2

w
*layerfilter1_neib_att_conv_head_1/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

)layerfilter1_neib_att_conv_head_1/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:

0layerfilter1_neib_att_conv_head_1/bn/beta/AssignAssign)layerfilter1_neib_att_conv_head_1/bn/beta*layerfilter1_neib_att_conv_head_1/bn/Const*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/bn/beta*
_output_shapes
:
È
.layerfilter1_neib_att_conv_head_1/bn/beta/readIdentity)layerfilter1_neib_att_conv_head_1/bn/beta*
T0*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/bn/beta*
_output_shapes
:
y
,layerfilter1_neib_att_conv_head_1/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

*layerfilter1_neib_att_conv_head_1/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
¢
1layerfilter1_neib_att_conv_head_1/bn/gamma/AssignAssign*layerfilter1_neib_att_conv_head_1/bn/gamma,layerfilter1_neib_att_conv_head_1/bn/Const_1*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_1/bn/gamma*
_output_shapes
:
Ë
/layerfilter1_neib_att_conv_head_1/bn/gamma/readIdentity*layerfilter1_neib_att_conv_head_1/bn/gamma*
T0*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_1/bn/gamma*
_output_shapes
:

Clayerfilter1_neib_att_conv_head_1/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
÷
1layerfilter1_neib_att_conv_head_1/bn/moments/meanMean)layerfilter1_neib_att_conv_head_1/BiasAddClayerfilter1_neib_att_conv_head_1/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
­
9layerfilter1_neib_att_conv_head_1/bn/moments/StopGradientStopGradient1layerfilter1_neib_att_conv_head_1/bn/moments/mean*
T0*&
_output_shapes
:
ê
>layerfilter1_neib_att_conv_head_1/bn/moments/SquaredDifferenceSquaredDifference)layerfilter1_neib_att_conv_head_1/BiasAdd9layerfilter1_neib_att_conv_head_1/bn/moments/StopGradient*
T0*&
_output_shapes
:2


Glayerfilter1_neib_att_conv_head_1/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

5layerfilter1_neib_att_conv_head_1/bn/moments/varianceMean>layerfilter1_neib_att_conv_head_1/bn/moments/SquaredDifferenceGlayerfilter1_neib_att_conv_head_1/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:
°
4layerfilter1_neib_att_conv_head_1/bn/moments/SqueezeSqueeze1layerfilter1_neib_att_conv_head_1/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:
¶
6layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1Squeeze5layerfilter1_neib_att_conv_head_1/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:
{
0layerfilter1_neib_att_conv_head_1/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

2layerfilter1_neib_att_conv_head_1/bn/cond/switch_tIdentity2layerfilter1_neib_att_conv_head_1/bn/cond/Switch:1*
T0
*
_output_shapes
: 

2layerfilter1_neib_att_conv_head_1/bn/cond/switch_fIdentity0layerfilter1_neib_att_conv_head_1/bn/cond/Switch*
T0
*
_output_shapes
: 
m
1layerfilter1_neib_att_conv_head_1/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
ç
layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ø
layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ð
layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
æ
rlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
ylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragelayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¤
wlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ë
layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ø
layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilllayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ê
tlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
{layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssigntlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ª
ylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentitytlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Â
Hlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/decayConst3^layerfilter1_neib_att_conv_head_1/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
Ú
Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst3^layerfilter1_neib_att_conv_head_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ª
Vlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubXlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xHlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/decay*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ô
Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Subalayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1clayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
È
_layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchwlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read1layerfilter1_neib_att_conv_head_1/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
È
alayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch4layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze1layerfilter1_neib_att_conv_head_1/bn/cond/pred_id*
T0*G
_class=
;9loc:@layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze* 
_output_shapes
::
¼
Vlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulXlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Vlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ô
Rlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub[layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Vlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
À
Ylayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage1layerfilter1_neib_att_conv_head_1/bn/cond/pred_id*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
Þ
Zlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst3^layerfilter1_neib_att_conv_head_1/bn/cond/switch_t*
valueB
 *  ?*
dtype0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
°
Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubZlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xHlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/decay*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ü
Zlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Subclayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1elayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Î
alayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read1layerfilter1_neib_att_conv_head_1/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Î
clayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch6layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_11layerfilter1_neib_att_conv_head_1/bn/cond/pred_id*
T0*I
_class?
=;loc:@layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1* 
_output_shapes
::
Ä
Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulZlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ü
Tlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub]layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Æ
[layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchtlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage1layerfilter1_neib_att_conv_head_1/bn/cond/pred_id*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
«
Blayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverageNoOpS^layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvgU^layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_13^layerfilter1_neib_att_conv_head_1/bn/cond/switch_t
©
<layerfilter1_neib_att_conv_head_1/bn/cond/control_dependencyIdentity2layerfilter1_neib_att_conv_head_1/bn/cond/switch_tC^layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage*
T0
*E
_class;
97loc:@layerfilter1_neib_att_conv_head_1/bn/cond/switch_t*
_output_shapes
: 
k
.layerfilter1_neib_att_conv_head_1/bn/cond/NoOpNoOp3^layerfilter1_neib_att_conv_head_1/bn/cond/switch_f

>layerfilter1_neib_att_conv_head_1/bn/cond/control_dependency_1Identity2layerfilter1_neib_att_conv_head_1/bn/cond/switch_f/^layerfilter1_neib_att_conv_head_1/bn/cond/NoOp*
T0
*E
_class;
97loc:@layerfilter1_neib_att_conv_head_1/bn/cond/switch_f*
_output_shapes
: 
â
/layerfilter1_neib_att_conv_head_1/bn/cond/MergeMerge>layerfilter1_neib_att_conv_head_1/bn/cond/control_dependency_1<layerfilter1_neib_att_conv_head_1/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
}
2layerfilter1_neib_att_conv_head_1/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 

4layerfilter1_neib_att_conv_head_1/bn/cond_1/switch_tIdentity4layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 

4layerfilter1_neib_att_conv_head_1/bn/cond_1/switch_fIdentity2layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch*
T0
*
_output_shapes
: 
o
3layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
Ö
4layerfilter1_neib_att_conv_head_1/bn/cond_1/IdentityIdentity=layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity/Switch:10^layerfilter1_neib_att_conv_head_1/bn/cond/Merge*
T0*
_output_shapes
:
¤
;layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity/SwitchSwitch4layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze3layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id*
T0*G
_class=
;9loc:@layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze* 
_output_shapes
::
Ú
6layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity_1Identity?layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity_1/Switch:10^layerfilter1_neib_att_conv_head_1/bn/cond/Merge*
T0*
_output_shapes
:
ª
=layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity_1/SwitchSwitch6layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_13layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1* 
_output_shapes
::

4layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_1Switchwlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read3layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id*
T0*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
£
4layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_2Switchylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read3layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id*
T0*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Ö
1layerfilter1_neib_att_conv_head_1/bn/cond_1/MergeMerge4layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_14layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 
Ú
3layerfilter1_neib_att_conv_head_1/bn/cond_1/Merge_1Merge4layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_26layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
y
4layerfilter1_neib_att_conv_head_1/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
É
2layerfilter1_neib_att_conv_head_1/bn/batchnorm/addAdd3layerfilter1_neib_att_conv_head_1/bn/cond_1/Merge_14layerfilter1_neib_att_conv_head_1/bn/batchnorm/add/y*
T0*
_output_shapes
:

4layerfilter1_neib_att_conv_head_1/bn/batchnorm/RsqrtRsqrt2layerfilter1_neib_att_conv_head_1/bn/batchnorm/add*
T0*
_output_shapes
:
Å
2layerfilter1_neib_att_conv_head_1/bn/batchnorm/mulMul4layerfilter1_neib_att_conv_head_1/bn/batchnorm/Rsqrt/layerfilter1_neib_att_conv_head_1/bn/gamma/read*
T0*
_output_shapes
:
Ë
4layerfilter1_neib_att_conv_head_1/bn/batchnorm/mul_1Mul)layerfilter1_neib_att_conv_head_1/BiasAdd2layerfilter1_neib_att_conv_head_1/bn/batchnorm/mul*
T0*&
_output_shapes
:2

Ç
4layerfilter1_neib_att_conv_head_1/bn/batchnorm/mul_2Mul1layerfilter1_neib_att_conv_head_1/bn/cond_1/Merge2layerfilter1_neib_att_conv_head_1/bn/batchnorm/mul*
T0*
_output_shapes
:
Ä
2layerfilter1_neib_att_conv_head_1/bn/batchnorm/subSub.layerfilter1_neib_att_conv_head_1/bn/beta/read4layerfilter1_neib_att_conv_head_1/bn/batchnorm/mul_2*
T0*
_output_shapes
:
Ö
4layerfilter1_neib_att_conv_head_1/bn/batchnorm/add_1Add4layerfilter1_neib_att_conv_head_1/bn/batchnorm/mul_12layerfilter1_neib_att_conv_head_1/bn/batchnorm/sub*
T0*&
_output_shapes
:2


&layerfilter1_neib_att_conv_head_1/ReluRelu4layerfilter1_neib_att_conv_head_1/bn/batchnorm/add_1*
T0*&
_output_shapes
:2


add_12Add&layerfilter1_self_att_conv_head_1/Relu&layerfilter1_neib_att_conv_head_1/Relu*
T0*&
_output_shapes
:2

i
transpose_8/permConst*%
valueB"             *
dtype0*
_output_shapes
:
p
transpose_8	Transposeadd_12transpose_8/perm*
T0*
Tperm0*&
_output_shapes
:2

V
LeakyRelu_2/alphaConst*
valueB
 *ÍÌL>*
dtype0*
_output_shapes
: 
g
LeakyRelu_2/mulMulLeakyRelu_2/alphatranspose_8*
T0*&
_output_shapes
:2

m
LeakyRelu_2/MaximumMaximumLeakyRelu_2/multranspose_8*
T0*&
_output_shapes
:2

`
Shape_5Const*%
valueB"   2      
   *
dtype0*
_output_shapes
:
H
Rank_2Const*
value	B :*
dtype0*
_output_shapes
: 
`
Shape_6Const*%
valueB"   2      
   *
dtype0*
_output_shapes
:
I
Sub_2/yConst*
value	B :*
dtype0*
_output_shapes
: 
>
Sub_2SubRank_2Sub_2/y*
T0*
_output_shapes
: 
V
Slice_2/beginPackSub_2*
N*
T0*

axis *
_output_shapes
:
V
Slice_2/sizeConst*
valueB:*
dtype0*
_output_shapes
:
h
Slice_2SliceShape_6Slice_2/beginSlice_2/size*
T0*
Index0*
_output_shapes
:
d
concat_5/values_0Const*
valueB:
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
:
O
concat_5/axisConst*
value	B : *
dtype0*
_output_shapes
: 
y
concat_5ConcatV2concat_5/values_0Slice_2concat_5/axis*
N*
T0*

Tidx0*
_output_shapes
:
k

Reshape_10ReshapeLeakyRelu_2/Maximumconcat_5*
T0*
Tshape0*
_output_shapes

:2

I
	Softmax_2Softmax
Reshape_10*
T0*
_output_shapes

:2

h

Reshape_11Reshape	Softmax_2Shape_5*
T0*
Tshape0*&
_output_shapes
:2


MatMul_4BatchMatMul
Reshape_11layerfilter1_edgefea_1/Relu*
T0*
adj_x( *
adj_y( *&
_output_shapes
:2@
¡
2BiasAdd_2/biases/Initializer/zeros/shape_as_tensorConst*
valueB:@*
dtype0*#
_class
loc:@BiasAdd_2/biases*
_output_shapes
:

(BiasAdd_2/biases/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*#
_class
loc:@BiasAdd_2/biases*
_output_shapes
: 
ä
"BiasAdd_2/biases/Initializer/zerosFill2BiasAdd_2/biases/Initializer/zeros/shape_as_tensor(BiasAdd_2/biases/Initializer/zeros/Const*
T0*

index_type0*#
_class
loc:@BiasAdd_2/biases*
_output_shapes
:@
¡
BiasAdd_2/biases
VariableV2*
shape:@*
dtype0*
	container *
shared_name *#
_class
loc:@BiasAdd_2/biases*
_output_shapes
:@
Ê
BiasAdd_2/biases/AssignAssignBiasAdd_2/biases"BiasAdd_2/biases/Initializer/zeros*
T0*
validate_shape(*
use_locking(*#
_class
loc:@BiasAdd_2/biases*
_output_shapes
:@
}
BiasAdd_2/biases/readIdentityBiasAdd_2/biases*
T0*#
_class
loc:@BiasAdd_2/biases*
_output_shapes
:@

BiasAdd_2/BiasAddBiasAddMatMul_4BiasAdd_2/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2@
R
Relu_2ReluBiasAdd_2/BiasAdd*
T0*&
_output_shapes
:2@
X
concat_6/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
z
concat_6ConcatV2Relu_1Relu_2concat_6/axis*
N*
T0*

Tidx0*'
_output_shapes
:2
\
ExpandDims_15/dimConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
x
ExpandDims_15
ExpandDimsPlaceholderExpandDims_15/dim*
T0*

Tdim0*&
_output_shapes
:2
X
concat_7/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 

concat_7ConcatV2ExpandDims_15concat_6concat_7/axis*
N*
T0*

Tidx0*'
_output_shapes
:2
X
concat_8/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
¤
concat_8ConcatV2layerfilter1_edgefea_0/Relulayerfilter1_edgefea_1/Reluconcat_8/axis*
N*
T0*

Tidx0*'
_output_shapes
:2

b
Max_1/reduction_indicesConst*
valueB :
þÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
~
Max_1Maxconcat_8Max_1/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:2
¯
1gapnet10/weights/Initializer/random_uniform/shapeConst*%
valueB"            *
dtype0*#
_class
loc:@gapnet10/weights*
_output_shapes
:

/gapnet10/weights/Initializer/random_uniform/minConst*
valueB
 *Ñ_ý½*
dtype0*#
_class
loc:@gapnet10/weights*
_output_shapes
: 

/gapnet10/weights/Initializer/random_uniform/maxConst*
valueB
 *Ñ_ý=*
dtype0*#
_class
loc:@gapnet10/weights*
_output_shapes
: 
û
9gapnet10/weights/Initializer/random_uniform/RandomUniformRandomUniform1gapnet10/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*#
_class
loc:@gapnet10/weights*(
_output_shapes
:
Þ
/gapnet10/weights/Initializer/random_uniform/subSub/gapnet10/weights/Initializer/random_uniform/max/gapnet10/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet10/weights*
_output_shapes
: 
ú
/gapnet10/weights/Initializer/random_uniform/mulMul9gapnet10/weights/Initializer/random_uniform/RandomUniform/gapnet10/weights/Initializer/random_uniform/sub*
T0*#
_class
loc:@gapnet10/weights*(
_output_shapes
:
ì
+gapnet10/weights/Initializer/random_uniformAdd/gapnet10/weights/Initializer/random_uniform/mul/gapnet10/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet10/weights*(
_output_shapes
:
Ì
gapnet10/weights
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *#
_class
loc:@gapnet10/weights*(
_output_shapes
:
ð
gapnet10/weights/AssignAssigngapnet10/weights+gapnet10/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet10/weights*(
_output_shapes
:

gapnet10/weights/readIdentitygapnet10/weights"/device:CPU:0*
T0*#
_class
loc:@gapnet10/weights*(
_output_shapes
:
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
Û
gapnet10/Conv2DConv2Dconcat_7gapnet10/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*'
_output_shapes
:2

!gapnet10/biases/Initializer/ConstConst*
valueB*    *
dtype0*"
_class
loc:@gapnet10/biases*
_output_shapes	
:
°
gapnet10/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *"
_class
loc:@gapnet10/biases*
_output_shapes	
:
Ö
gapnet10/biases/AssignAssigngapnet10/biases!gapnet10/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet10/biases*
_output_shapes	
:

gapnet10/biases/readIdentitygapnet10/biases"/device:CPU:0*
T0*"
_class
loc:@gapnet10/biases*
_output_shapes	
:

gapnet10/BiasAddBiasAddgapnet10/Conv2Dgapnet10/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:2
`
gapnet10/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
~
gapnet10/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
º
gapnet10/bn/beta/AssignAssigngapnet10/bn/betagapnet10/bn/Const*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet10/bn/beta*
_output_shapes	
:
~
gapnet10/bn/beta/readIdentitygapnet10/bn/beta*
T0*#
_class
loc:@gapnet10/bn/beta*
_output_shapes	
:
b
gapnet10/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:

gapnet10/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
¿
gapnet10/bn/gamma/AssignAssigngapnet10/bn/gammagapnet10/bn/Const_1*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet10/bn/gamma*
_output_shapes	
:

gapnet10/bn/gamma/readIdentitygapnet10/bn/gamma*
T0*$
_class
loc:@gapnet10/bn/gamma*
_output_shapes	
:

*gapnet10/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
­
gapnet10/bn/moments/meanMeangapnet10/BiasAdd*gapnet10/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:
|
 gapnet10/bn/moments/StopGradientStopGradientgapnet10/bn/moments/mean*
T0*'
_output_shapes
:
 
%gapnet10/bn/moments/SquaredDifferenceSquaredDifferencegapnet10/BiasAdd gapnet10/bn/moments/StopGradient*
T0*'
_output_shapes
:2

.gapnet10/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ê
gapnet10/bn/moments/varianceMean%gapnet10/bn/moments/SquaredDifference.gapnet10/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:

gapnet10/bn/moments/SqueezeSqueezegapnet10/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes	
:

gapnet10/bn/moments/Squeeze_1Squeezegapnet10/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes	
:
b
gapnet10/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
a
gapnet10/bn/cond/switch_tIdentitygapnet10/bn/cond/Switch:1*
T0
*
_output_shapes
: 
_
gapnet10/bn/cond/switch_fIdentitygapnet10/bn/cond/Switch*
T0
*
_output_shapes
: 
T
gapnet10/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

bgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ò
Xgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¥
Rgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFillbgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorXgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageRgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Egapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

dgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ö
Zgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
­
Tgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilldgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorZgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Bgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageTgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Ggapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

/gapnet10/bn/cond/ExponentialMovingAverage/decayConst^gapnet10/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
õ
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet10/bn/cond/switch_t*
valueB
 *  ?*
dtype0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¬
=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x/gapnet10/bn/cond/ExponentialMovingAverage/decay*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
×
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubHgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
³
Fgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchEgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet10/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
æ
Hgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet10/bn/moments/Squeezegapnet10/bn/cond/pred_id*
T0*.
_class$
" loc:@gapnet10/bn/moments/Squeeze*"
_output_shapes
::
¿
=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
×
9gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubBgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1=gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
«
@gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAveragegapnet10/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
ù
Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet10/bn/cond/switch_t*
valueB
 *  ?*
dtype0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
²
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubAgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x/gapnet10/bn/cond/ExponentialMovingAverage/decay*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ß
Agapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubJgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Lgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
¹
Hgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchGgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet10/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
ì
Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchgapnet10/bn/moments/Squeeze_1gapnet10/bn/cond/pred_id*
T0*0
_class&
$"loc:@gapnet10/bn/moments/Squeeze_1*"
_output_shapes
::
Ç
?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulAgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
ß
;gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubDgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1?gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
±
Bgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet10/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
Ç
)gapnet10/bn/cond/ExponentialMovingAverageNoOp:^gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg<^gapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1^gapnet10/bn/cond/switch_t
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
%gapnet10/bn/cond/control_dependency_1Identitygapnet10/bn/cond/switch_f^gapnet10/bn/cond/NoOp*
T0
*,
_class"
 loc:@gapnet10/bn/cond/switch_f*
_output_shapes
: 

gapnet10/bn/cond/MergeMerge%gapnet10/bn/cond/control_dependency_1#gapnet10/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
d
gapnet10/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
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
gapnet10/bn/cond_1/switch_fIdentitygapnet10/bn/cond_1/Switch*
T0
*
_output_shapes
: 
V
gapnet10/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

gapnet10/bn/cond_1/IdentityIdentity$gapnet10/bn/cond_1/Identity/Switch:1^gapnet10/bn/cond/Merge*
T0*
_output_shapes	
:
Â
"gapnet10/bn/cond_1/Identity/SwitchSwitchgapnet10/bn/moments/Squeezegapnet10/bn/cond_1/pred_id*
T0*.
_class$
" loc:@gapnet10/bn/moments/Squeeze*"
_output_shapes
::

gapnet10/bn/cond_1/Identity_1Identity&gapnet10/bn/cond_1/Identity_1/Switch:1^gapnet10/bn/cond/Merge*
T0*
_output_shapes	
:
È
$gapnet10/bn/cond_1/Identity_1/SwitchSwitchgapnet10/bn/moments/Squeeze_1gapnet10/bn/cond_1/pred_id*
T0*0
_class&
$"loc:@gapnet10/bn/moments/Squeeze_1*"
_output_shapes
::

gapnet10/bn/cond_1/Switch_1SwitchEgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet10/bn/cond_1/pred_id*
T0*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

gapnet10/bn/cond_1/Switch_2SwitchGgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet10/bn/cond_1/pred_id*
T0*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::

gapnet10/bn/cond_1/MergeMergegapnet10/bn/cond_1/Switch_1gapnet10/bn/cond_1/Identity*
T0*
N*
_output_shapes
	:: 

gapnet10/bn/cond_1/Merge_1Mergegapnet10/bn/cond_1/Switch_2gapnet10/bn/cond_1/Identity_1*
T0*
N*
_output_shapes
	:: 
`
gapnet10/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 

gapnet10/bn/batchnorm/addAddgapnet10/bn/cond_1/Merge_1gapnet10/bn/batchnorm/add/y*
T0*
_output_shapes	
:
e
gapnet10/bn/batchnorm/RsqrtRsqrtgapnet10/bn/batchnorm/add*
T0*
_output_shapes	
:
{
gapnet10/bn/batchnorm/mulMulgapnet10/bn/batchnorm/Rsqrtgapnet10/bn/gamma/read*
T0*
_output_shapes	
:

gapnet10/bn/batchnorm/mul_1Mulgapnet10/BiasAddgapnet10/bn/batchnorm/mul*
T0*'
_output_shapes
:2
}
gapnet10/bn/batchnorm/mul_2Mulgapnet10/bn/cond_1/Mergegapnet10/bn/batchnorm/mul*
T0*
_output_shapes	
:
z
gapnet10/bn/batchnorm/subSubgapnet10/bn/beta/readgapnet10/bn/batchnorm/mul_2*
T0*
_output_shapes	
:

gapnet10/bn/batchnorm/add_1Addgapnet10/bn/batchnorm/mul_1gapnet10/bn/batchnorm/sub*
T0*'
_output_shapes
:2
d
gapnet10/ReluRelugapnet10/bn/batchnorm/add_1*
T0*'
_output_shapes
:2
¯
1gapnet11/weights/Initializer/random_uniform/shapeConst*%
valueB"            *
dtype0*#
_class
loc:@gapnet11/weights*
_output_shapes
:

/gapnet11/weights/Initializer/random_uniform/minConst*
valueB
 *   ¾*
dtype0*#
_class
loc:@gapnet11/weights*
_output_shapes
: 

/gapnet11/weights/Initializer/random_uniform/maxConst*
valueB
 *   >*
dtype0*#
_class
loc:@gapnet11/weights*
_output_shapes
: 
û
9gapnet11/weights/Initializer/random_uniform/RandomUniformRandomUniform1gapnet11/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*#
_class
loc:@gapnet11/weights*(
_output_shapes
:
Þ
/gapnet11/weights/Initializer/random_uniform/subSub/gapnet11/weights/Initializer/random_uniform/max/gapnet11/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet11/weights*
_output_shapes
: 
ú
/gapnet11/weights/Initializer/random_uniform/mulMul9gapnet11/weights/Initializer/random_uniform/RandomUniform/gapnet11/weights/Initializer/random_uniform/sub*
T0*#
_class
loc:@gapnet11/weights*(
_output_shapes
:
ì
+gapnet11/weights/Initializer/random_uniformAdd/gapnet11/weights/Initializer/random_uniform/mul/gapnet11/weights/Initializer/random_uniform/min*
T0*#
_class
loc:@gapnet11/weights*(
_output_shapes
:
Ì
gapnet11/weights
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *#
_class
loc:@gapnet11/weights*(
_output_shapes
:
ð
gapnet11/weights/AssignAssigngapnet11/weights+gapnet11/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet11/weights*(
_output_shapes
:

gapnet11/weights/readIdentitygapnet11/weights"/device:CPU:0*
T0*#
_class
loc:@gapnet11/weights*(
_output_shapes
:
Q
gapnet11/L2LossL2Lossgapnet11/weights/read*
T0*
_output_shapes
: 
[
gapnet11/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
e
gapnet11/weight_lossMulgapnet11/L2Lossgapnet11/weight_loss/y*
T0*
_output_shapes
: 
à
gapnet11/Conv2DConv2Dgapnet10/Relugapnet11/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*'
_output_shapes
:2

!gapnet11/biases/Initializer/ConstConst*
valueB*    *
dtype0*"
_class
loc:@gapnet11/biases*
_output_shapes	
:
°
gapnet11/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *"
_class
loc:@gapnet11/biases*
_output_shapes	
:
Ö
gapnet11/biases/AssignAssigngapnet11/biases!gapnet11/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet11/biases*
_output_shapes	
:

gapnet11/biases/readIdentitygapnet11/biases"/device:CPU:0*
T0*"
_class
loc:@gapnet11/biases*
_output_shapes	
:

gapnet11/BiasAddBiasAddgapnet11/Conv2Dgapnet11/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:2
`
gapnet11/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
~
gapnet11/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
º
gapnet11/bn/beta/AssignAssigngapnet11/bn/betagapnet11/bn/Const*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet11/bn/beta*
_output_shapes	
:
~
gapnet11/bn/beta/readIdentitygapnet11/bn/beta*
T0*#
_class
loc:@gapnet11/bn/beta*
_output_shapes	
:
b
gapnet11/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:

gapnet11/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
¿
gapnet11/bn/gamma/AssignAssigngapnet11/bn/gammagapnet11/bn/Const_1*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet11/bn/gamma*
_output_shapes	
:

gapnet11/bn/gamma/readIdentitygapnet11/bn/gamma*
T0*$
_class
loc:@gapnet11/bn/gamma*
_output_shapes	
:

*gapnet11/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
­
gapnet11/bn/moments/meanMeangapnet11/BiasAdd*gapnet11/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:
|
 gapnet11/bn/moments/StopGradientStopGradientgapnet11/bn/moments/mean*
T0*'
_output_shapes
:
 
%gapnet11/bn/moments/SquaredDifferenceSquaredDifferencegapnet11/BiasAdd gapnet11/bn/moments/StopGradient*
T0*'
_output_shapes
:2

.gapnet11/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ê
gapnet11/bn/moments/varianceMean%gapnet11/bn/moments/SquaredDifference.gapnet11/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:

gapnet11/bn/moments/SqueezeSqueezegapnet11/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes	
:

gapnet11/bn/moments/Squeeze_1Squeezegapnet11/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes	
:
b
gapnet11/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
a
gapnet11/bn/cond/switch_tIdentitygapnet11/bn/cond/Switch:1*
T0
*
_output_shapes
: 
_
gapnet11/bn/cond/switch_fIdentitygapnet11/bn/cond/Switch*
T0
*
_output_shapes
: 
T
gapnet11/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

bgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ò
Xgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¥
Rgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFillbgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorXgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Ggapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverageRgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Egapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
T0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

dgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ö
Zgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
­
Tgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFilldgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorZgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Bgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Igapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverageTgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Ggapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

/gapnet11/bn/cond/ExponentialMovingAverage/decayConst^gapnet11/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
õ
?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^gapnet11/bn/cond/switch_t*
valueB
 *  ?*
dtype0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¬
=gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x/gapnet11/bn/cond/ExponentialMovingAverage/decay*
T0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
×
?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubHgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Jgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
³
Fgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchEgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet11/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
æ
Hgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchgapnet11/bn/moments/Squeezegapnet11/bn/cond/pred_id*
T0*.
_class$
" loc:@gapnet11/bn/moments/Squeeze*"
_output_shapes
::
¿
=gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1=gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
×
9gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubBgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1=gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
«
@gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAveragegapnet11/bn/cond/pred_id*
T0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
ù
Agapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^gapnet11/bn/cond/switch_t*
valueB
 *  ?*
dtype0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
²
?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubAgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x/gapnet11/bn/cond/ExponentialMovingAverage/decay*
T0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ß
Agapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubJgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Lgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
¹
Hgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchGgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet11/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
ì
Jgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchgapnet11/bn/moments/Squeeze_1gapnet11/bn/cond/pred_id*
T0*0
_class&
$"loc:@gapnet11/bn/moments/Squeeze_1*"
_output_shapes
::
Ç
?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulAgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
ß
;gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubDgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
±
Bgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet11/bn/cond/pred_id*
T0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
Ç
)gapnet11/bn/cond/ExponentialMovingAverageNoOp:^gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg<^gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1^gapnet11/bn/cond/switch_t
Å
#gapnet11/bn/cond/control_dependencyIdentitygapnet11/bn/cond/switch_t*^gapnet11/bn/cond/ExponentialMovingAverage*
T0
*,
_class"
 loc:@gapnet11/bn/cond/switch_t*
_output_shapes
: 
9
gapnet11/bn/cond/NoOpNoOp^gapnet11/bn/cond/switch_f
³
%gapnet11/bn/cond/control_dependency_1Identitygapnet11/bn/cond/switch_f^gapnet11/bn/cond/NoOp*
T0
*,
_class"
 loc:@gapnet11/bn/cond/switch_f*
_output_shapes
: 

gapnet11/bn/cond/MergeMerge%gapnet11/bn/cond/control_dependency_1#gapnet11/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
d
gapnet11/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
e
gapnet11/bn/cond_1/switch_tIdentitygapnet11/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
c
gapnet11/bn/cond_1/switch_fIdentitygapnet11/bn/cond_1/Switch*
T0
*
_output_shapes
: 
V
gapnet11/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

gapnet11/bn/cond_1/IdentityIdentity$gapnet11/bn/cond_1/Identity/Switch:1^gapnet11/bn/cond/Merge*
T0*
_output_shapes	
:
Â
"gapnet11/bn/cond_1/Identity/SwitchSwitchgapnet11/bn/moments/Squeezegapnet11/bn/cond_1/pred_id*
T0*.
_class$
" loc:@gapnet11/bn/moments/Squeeze*"
_output_shapes
::

gapnet11/bn/cond_1/Identity_1Identity&gapnet11/bn/cond_1/Identity_1/Switch:1^gapnet11/bn/cond/Merge*
T0*
_output_shapes	
:
È
$gapnet11/bn/cond_1/Identity_1/SwitchSwitchgapnet11/bn/moments/Squeeze_1gapnet11/bn/cond_1/pred_id*
T0*0
_class&
$"loc:@gapnet11/bn/moments/Squeeze_1*"
_output_shapes
::

gapnet11/bn/cond_1/Switch_1SwitchEgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/readgapnet11/bn/cond_1/pred_id*
T0*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::

gapnet11/bn/cond_1/Switch_2SwitchGgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/readgapnet11/bn/cond_1/pred_id*
T0*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::

gapnet11/bn/cond_1/MergeMergegapnet11/bn/cond_1/Switch_1gapnet11/bn/cond_1/Identity*
T0*
N*
_output_shapes
	:: 

gapnet11/bn/cond_1/Merge_1Mergegapnet11/bn/cond_1/Switch_2gapnet11/bn/cond_1/Identity_1*
T0*
N*
_output_shapes
	:: 
`
gapnet11/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 

gapnet11/bn/batchnorm/addAddgapnet11/bn/cond_1/Merge_1gapnet11/bn/batchnorm/add/y*
T0*
_output_shapes	
:
e
gapnet11/bn/batchnorm/RsqrtRsqrtgapnet11/bn/batchnorm/add*
T0*
_output_shapes	
:
{
gapnet11/bn/batchnorm/mulMulgapnet11/bn/batchnorm/Rsqrtgapnet11/bn/gamma/read*
T0*
_output_shapes	
:

gapnet11/bn/batchnorm/mul_1Mulgapnet11/BiasAddgapnet11/bn/batchnorm/mul*
T0*'
_output_shapes
:2
}
gapnet11/bn/batchnorm/mul_2Mulgapnet11/bn/cond_1/Mergegapnet11/bn/batchnorm/mul*
T0*
_output_shapes	
:
z
gapnet11/bn/batchnorm/subSubgapnet11/bn/beta/readgapnet11/bn/batchnorm/mul_2*
T0*
_output_shapes	
:

gapnet11/bn/batchnorm/add_1Addgapnet11/bn/batchnorm/mul_1gapnet11/bn/batchnorm/sub*
T0*'
_output_shapes
:2
d
gapnet11/ReluRelugapnet11/bn/batchnorm/add_1*
T0*'
_output_shapes
:2
i
Reshape_12/shapeConst*%
valueB"         ÿÿÿÿ*
dtype0*
_output_shapes
:
u

Reshape_12ReshapePlaceholder_2Reshape_12/shape*
T0*
Tshape0*&
_output_shapes
:
i
Tile_4/multiplesConst*%
valueB"   2         *
dtype0*
_output_shapes
:
o
Tile_4Tile
Reshape_12Tile_4/multiples*
T0*

Tmultiples0*&
_output_shapes
:2
¹
6global_expand/weights/Initializer/random_uniform/shapeConst*%
valueB"            *
dtype0*(
_class
loc:@global_expand/weights*
_output_shapes
:
£
4global_expand/weights/Initializer/random_uniform/minConst*
valueB
 *:Í¿*
dtype0*(
_class
loc:@global_expand/weights*
_output_shapes
: 
£
4global_expand/weights/Initializer/random_uniform/maxConst*
valueB
 *:Í?*
dtype0*(
_class
loc:@global_expand/weights*
_output_shapes
: 

>global_expand/weights/Initializer/random_uniform/RandomUniformRandomUniform6global_expand/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*(
_class
loc:@global_expand/weights*&
_output_shapes
:
ò
4global_expand/weights/Initializer/random_uniform/subSub4global_expand/weights/Initializer/random_uniform/max4global_expand/weights/Initializer/random_uniform/min*
T0*(
_class
loc:@global_expand/weights*
_output_shapes
: 

4global_expand/weights/Initializer/random_uniform/mulMul>global_expand/weights/Initializer/random_uniform/RandomUniform4global_expand/weights/Initializer/random_uniform/sub*
T0*(
_class
loc:@global_expand/weights*&
_output_shapes
:
þ
0global_expand/weights/Initializer/random_uniformAdd4global_expand/weights/Initializer/random_uniform/mul4global_expand/weights/Initializer/random_uniform/min*
T0*(
_class
loc:@global_expand/weights*&
_output_shapes
:
Ò
global_expand/weights
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *(
_class
loc:@global_expand/weights*&
_output_shapes
:

global_expand/weights/AssignAssignglobal_expand/weights0global_expand/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@global_expand/weights*&
_output_shapes
:
§
global_expand/weights/readIdentityglobal_expand/weights"/device:CPU:0*
T0*(
_class
loc:@global_expand/weights*&
_output_shapes
:
[
global_expand/L2LossL2Lossglobal_expand/weights/read*
T0*
_output_shapes
: 
`
global_expand/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
t
global_expand/weight_lossMulglobal_expand/L2Lossglobal_expand/weight_loss/y*
T0*
_output_shapes
: 
â
global_expand/Conv2DConv2DTile_4global_expand/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2

&global_expand/biases/Initializer/ConstConst*
valueB*    *
dtype0*'
_class
loc:@global_expand/biases*
_output_shapes
:
¸
global_expand/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *'
_class
loc:@global_expand/biases*
_output_shapes
:
é
global_expand/biases/AssignAssignglobal_expand/biases&global_expand/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@global_expand/biases*
_output_shapes
:

global_expand/biases/readIdentityglobal_expand/biases"/device:CPU:0*
T0*'
_class
loc:@global_expand/biases*
_output_shapes
:

global_expand/BiasAddBiasAddglobal_expand/Conv2Dglobal_expand/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
c
global_expand/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes
:

global_expand/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
Í
global_expand/bn/beta/AssignAssignglobal_expand/bn/betaglobal_expand/bn/Const*
T0*
validate_shape(*
use_locking(*(
_class
loc:@global_expand/bn/beta*
_output_shapes
:

global_expand/bn/beta/readIdentityglobal_expand/bn/beta*
T0*(
_class
loc:@global_expand/bn/beta*
_output_shapes
:
e
global_expand/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes
:

global_expand/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes
:
Ò
global_expand/bn/gamma/AssignAssignglobal_expand/bn/gammaglobal_expand/bn/Const_1*
T0*
validate_shape(*
use_locking(*)
_class
loc:@global_expand/bn/gamma*
_output_shapes
:

global_expand/bn/gamma/readIdentityglobal_expand/bn/gamma*
T0*)
_class
loc:@global_expand/bn/gamma*
_output_shapes
:

/global_expand/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
»
global_expand/bn/moments/meanMeanglobal_expand/BiasAdd/global_expand/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:

%global_expand/bn/moments/StopGradientStopGradientglobal_expand/bn/moments/mean*
T0*&
_output_shapes
:
®
*global_expand/bn/moments/SquaredDifferenceSquaredDifferenceglobal_expand/BiasAdd%global_expand/bn/moments/StopGradient*
T0*&
_output_shapes
:2

3global_expand/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
Ø
!global_expand/bn/moments/varianceMean*global_expand/bn/moments/SquaredDifference3global_expand/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*&
_output_shapes
:

 global_expand/bn/moments/SqueezeSqueezeglobal_expand/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes
:

"global_expand/bn/moments/Squeeze_1Squeeze!global_expand/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes
:
g
global_expand/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
k
global_expand/bn/cond/switch_tIdentityglobal_expand/bn/cond/Switch:1*
T0
*
_output_shapes
: 
i
global_expand/bn/cond/switch_fIdentityglobal_expand/bn/cond/Switch*
T0
*
_output_shapes
: 
Y
global_expand/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

lglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:

bglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ì
\global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFilllglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorbglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:

Jglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
²
Qglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssignJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage\global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
«
Oglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/readIdentityJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
T0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:

nglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:

dglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ô
^global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFillnglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensordglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:

Lglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
º
Sglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssignLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage^global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
±
Qglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentityLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:

4global_expand/bn/cond/ExponentialMovingAverage/decayConst^global_expand/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 

Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^global_expand/bn/cond/switch_t*
valueB
 *  ?*
dtype0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Å
Bglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSubDglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x4global_expand/bn/cond/ExponentialMovingAverage/decay*
T0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ï
Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubMglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Oglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ï
Kglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitchOglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/readglobal_expand/bn/cond/pred_id*
T0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
ø
Mglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switch global_expand/bn/moments/Squeezeglobal_expand/bn/cond/pred_id*
T0*3
_class)
'%loc:@global_expand/bn/moments/Squeeze* 
_output_shapes
::
×
Bglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMulDglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1Bglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ï
>global_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSubGglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Bglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Ç
Eglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitchJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverageglobal_expand/bn/cond/pred_id*
T0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::

Fglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^global_expand/bn/cond/switch_t*
valueB
 *  ?*
dtype0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Ë
Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSubFglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x4global_expand/bn/cond/ExponentialMovingAverage/decay*
T0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
÷
Fglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubOglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Qglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Õ
Mglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitchQglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/readglobal_expand/bn/cond/pred_id*
T0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
þ
Oglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switch"global_expand/bn/moments/Squeeze_1global_expand/bn/cond/pred_id*
T0*5
_class+
)'loc:@global_expand/bn/moments/Squeeze_1* 
_output_shapes
::
ß
Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMulFglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
÷
@global_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSubIglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Í
Gglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitchLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverageglobal_expand/bn/cond/pred_id*
T0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::
Û
.global_expand/bn/cond/ExponentialMovingAverageNoOp?^global_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvgA^global_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1^global_expand/bn/cond/switch_t
Ù
(global_expand/bn/cond/control_dependencyIdentityglobal_expand/bn/cond/switch_t/^global_expand/bn/cond/ExponentialMovingAverage*
T0
*1
_class'
%#loc:@global_expand/bn/cond/switch_t*
_output_shapes
: 
C
global_expand/bn/cond/NoOpNoOp^global_expand/bn/cond/switch_f
Ç
*global_expand/bn/cond/control_dependency_1Identityglobal_expand/bn/cond/switch_f^global_expand/bn/cond/NoOp*
T0
*1
_class'
%#loc:@global_expand/bn/cond/switch_f*
_output_shapes
: 
¦
global_expand/bn/cond/MergeMerge*global_expand/bn/cond/control_dependency_1(global_expand/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
i
global_expand/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
o
 global_expand/bn/cond_1/switch_tIdentity global_expand/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
m
 global_expand/bn/cond_1/switch_fIdentityglobal_expand/bn/cond_1/Switch*
T0
*
_output_shapes
: 
[
global_expand/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 

 global_expand/bn/cond_1/IdentityIdentity)global_expand/bn/cond_1/Identity/Switch:1^global_expand/bn/cond/Merge*
T0*
_output_shapes
:
Ô
'global_expand/bn/cond_1/Identity/SwitchSwitch global_expand/bn/moments/Squeezeglobal_expand/bn/cond_1/pred_id*
T0*3
_class)
'%loc:@global_expand/bn/moments/Squeeze* 
_output_shapes
::

"global_expand/bn/cond_1/Identity_1Identity+global_expand/bn/cond_1/Identity_1/Switch:1^global_expand/bn/cond/Merge*
T0*
_output_shapes
:
Ú
)global_expand/bn/cond_1/Identity_1/SwitchSwitch"global_expand/bn/moments/Squeeze_1global_expand/bn/cond_1/pred_id*
T0*5
_class+
)'loc:@global_expand/bn/moments/Squeeze_1* 
_output_shapes
::
¦
 global_expand/bn/cond_1/Switch_1SwitchOglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/readglobal_expand/bn/cond_1/pred_id*
T0*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage* 
_output_shapes
::
ª
 global_expand/bn/cond_1/Switch_2SwitchQglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/readglobal_expand/bn/cond_1/pred_id*
T0*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage* 
_output_shapes
::

global_expand/bn/cond_1/MergeMerge global_expand/bn/cond_1/Switch_1 global_expand/bn/cond_1/Identity*
T0*
N*
_output_shapes

:: 

global_expand/bn/cond_1/Merge_1Merge global_expand/bn/cond_1/Switch_2"global_expand/bn/cond_1/Identity_1*
T0*
N*
_output_shapes

:: 
e
 global_expand/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 

global_expand/bn/batchnorm/addAddglobal_expand/bn/cond_1/Merge_1 global_expand/bn/batchnorm/add/y*
T0*
_output_shapes
:
n
 global_expand/bn/batchnorm/RsqrtRsqrtglobal_expand/bn/batchnorm/add*
T0*
_output_shapes
:

global_expand/bn/batchnorm/mulMul global_expand/bn/batchnorm/Rsqrtglobal_expand/bn/gamma/read*
T0*
_output_shapes
:

 global_expand/bn/batchnorm/mul_1Mulglobal_expand/BiasAddglobal_expand/bn/batchnorm/mul*
T0*&
_output_shapes
:2

 global_expand/bn/batchnorm/mul_2Mulglobal_expand/bn/cond_1/Mergeglobal_expand/bn/batchnorm/mul*
T0*
_output_shapes
:

global_expand/bn/batchnorm/subSubglobal_expand/bn/beta/read global_expand/bn/batchnorm/mul_2*
T0*
_output_shapes
:

 global_expand/bn/batchnorm/add_1Add global_expand/bn/batchnorm/mul_1global_expand/bn/batchnorm/sub*
T0*&
_output_shapes
:2
m
global_expand/ReluRelu global_expand/bn/batchnorm/add_1*
T0*&
_output_shapes
:2
X
concat_9/axisConst*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
·
concat_9ConcatV2gapnet00/Relugapnet01/Relugapnet11/Reluglobal_expand/ReluMaxMax_1concat_9/axis*
N*
T0*

Tidx0*'
_output_shapes
:2ð
¥
,agg/weights/Initializer/random_uniform/shapeConst*%
valueB"      ð     *
dtype0*
_class
loc:@agg/weights*
_output_shapes
:

*agg/weights/Initializer/random_uniform/minConst*
valueB
 *Xï¶½*
dtype0*
_class
loc:@agg/weights*
_output_shapes
: 

*agg/weights/Initializer/random_uniform/maxConst*
valueB
 *Xï¶=*
dtype0*
_class
loc:@agg/weights*
_output_shapes
: 
ì
4agg/weights/Initializer/random_uniform/RandomUniformRandomUniform,agg/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*
_class
loc:@agg/weights*(
_output_shapes
:ð
Ê
*agg/weights/Initializer/random_uniform/subSub*agg/weights/Initializer/random_uniform/max*agg/weights/Initializer/random_uniform/min*
T0*
_class
loc:@agg/weights*
_output_shapes
: 
æ
*agg/weights/Initializer/random_uniform/mulMul4agg/weights/Initializer/random_uniform/RandomUniform*agg/weights/Initializer/random_uniform/sub*
T0*
_class
loc:@agg/weights*(
_output_shapes
:ð
Ø
&agg/weights/Initializer/random_uniformAdd*agg/weights/Initializer/random_uniform/mul*agg/weights/Initializer/random_uniform/min*
T0*
_class
loc:@agg/weights*(
_output_shapes
:ð
Â
agg/weights
VariableV2"/device:CPU:0*
shape:ð*
dtype0*
	container *
shared_name *
_class
loc:@agg/weights*(
_output_shapes
:ð
Ü
agg/weights/AssignAssignagg/weights&agg/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/weights*(
_output_shapes
:ð

agg/weights/readIdentityagg/weights"/device:CPU:0*
T0*
_class
loc:@agg/weights*(
_output_shapes
:ð
G

agg/L2LossL2Lossagg/weights/read*
T0*
_output_shapes
: 
V
agg/weight_loss/yConst*
valueB
 *    *
dtype0*
_output_shapes
: 
V
agg/weight_lossMul
agg/L2Lossagg/weight_loss/y*
T0*
_output_shapes
: 
Ñ

agg/Conv2DConv2Dconcat_9agg/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*'
_output_shapes
:2

agg/biases/Initializer/ConstConst*
valueB*    *
dtype0*
_class
loc:@agg/biases*
_output_shapes	
:
¦

agg/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *
_class
loc:@agg/biases*
_output_shapes	
:
Â
agg/biases/AssignAssign
agg/biasesagg/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/biases*
_output_shapes	
:
{
agg/biases/readIdentity
agg/biases"/device:CPU:0*
T0*
_class
loc:@agg/biases*
_output_shapes	
:
|
agg/BiasAddBiasAdd
agg/Conv2Dagg/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:2
[
agg/bn/ConstConst*
valueB*    *
dtype0*
_output_shapes	
:
y
agg/bn/beta
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
¦
agg/bn/beta/AssignAssignagg/bn/betaagg/bn/Const*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/bn/beta*
_output_shapes	
:
o
agg/bn/beta/readIdentityagg/bn/beta*
T0*
_class
loc:@agg/bn/beta*
_output_shapes	
:
]
agg/bn/Const_1Const*
valueB*  ?*
dtype0*
_output_shapes	
:
z
agg/bn/gamma
VariableV2*
shape:*
dtype0*
	container *
shared_name *
_output_shapes	
:
«
agg/bn/gamma/AssignAssignagg/bn/gammaagg/bn/Const_1*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/bn/gamma*
_output_shapes	
:
r
agg/bn/gamma/readIdentityagg/bn/gamma*
T0*
_class
loc:@agg/bn/gamma*
_output_shapes	
:
z
%agg/bn/moments/mean/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:

agg/bn/moments/meanMeanagg/BiasAdd%agg/bn/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:
r
agg/bn/moments/StopGradientStopGradientagg/bn/moments/mean*
T0*'
_output_shapes
:

 agg/bn/moments/SquaredDifferenceSquaredDifferenceagg/BiasAddagg/bn/moments/StopGradient*
T0*'
_output_shapes
:2
~
)agg/bn/moments/variance/reduction_indicesConst*!
valueB"          *
dtype0*
_output_shapes
:
»
agg/bn/moments/varianceMean agg/bn/moments/SquaredDifference)agg/bn/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:
u
agg/bn/moments/SqueezeSqueezeagg/bn/moments/mean*
T0*
squeeze_dims
 *
_output_shapes	
:
{
agg/bn/moments/Squeeze_1Squeezeagg/bn/moments/variance*
T0*
squeeze_dims
 *
_output_shapes	
:
]
agg/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
W
agg/bn/cond/switch_tIdentityagg/bn/cond/Switch:1*
T0
*
_output_shapes
: 
U
agg/bn/cond/switch_fIdentityagg/bn/cond/Switch*
T0
*
_output_shapes
: 
O
agg/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
î
Xagg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
Þ
Nagg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ý
Hagg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zerosFillXagg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorNagg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
ï
6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
ã
=agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/AssignAssign6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverageHagg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
ð
;agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/readIdentity6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
T0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
ò
Zagg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
â
Pagg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

Jagg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zerosFillZagg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/shape_as_tensorPagg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros/Const*
T0*

index_type0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
ó
8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage
VariableV2*
shape:*
dtype0*
	container *
shared_name *K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
ë
?agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignAssign8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageJagg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros*
T0*
validate_shape(*
use_locking(*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
ö
=agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/readIdentity8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
T0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

*agg/bn/cond/ExponentialMovingAverage/decayConst^agg/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
á
:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/xConst^agg/bn/cond/switch_t*
valueB
 *  ?*
dtype0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 

8agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/subSub:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x*agg/bn/cond/ExponentialMovingAverage/decay*
T0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
¾
:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1SubCagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1Eagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1*
T0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

Aagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/SwitchSwitch;agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/readagg/bn/cond/pred_id*
T0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
Ò
Cagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1Switchagg/bn/moments/Squeezeagg/bn/cond/pred_id*
T0*)
_class
loc:@agg/bn/moments/Squeeze*"
_output_shapes
::
¦
8agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mulMul:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_18agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub*
T0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
¾
4agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg	AssignSub=agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:18agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul*
T0*
use_locking( *I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

;agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch	RefSwitch6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverageagg/bn/cond/pred_id*
T0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
å
<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/xConst^agg/bn/cond/switch_t*
valueB
 *  ?*
dtype0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 

:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/subSub<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x*agg/bn/cond/ExponentialMovingAverage/decay*
T0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
Æ
<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1SubEagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1Gagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1*
T0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

Cagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/SwitchSwitch=agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/readagg/bn/cond/pred_id*
T0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
Ø
Eagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1Switchagg/bn/moments/Squeeze_1agg/bn/cond/pred_id*
T0*+
_class!
loc:@agg/bn/moments/Squeeze_1*"
_output_shapes
::
®
:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mulMul<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub*
T0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
Æ
6agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1	AssignSub?agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul*
T0*
use_locking( *K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:

=agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch	RefSwitch8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageagg/bn/cond/pred_id*
T0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
³
$agg/bn/cond/ExponentialMovingAverageNoOp5^agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg7^agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1^agg/bn/cond/switch_t
±
agg/bn/cond/control_dependencyIdentityagg/bn/cond/switch_t%^agg/bn/cond/ExponentialMovingAverage*
T0
*'
_class
loc:@agg/bn/cond/switch_t*
_output_shapes
: 
/
agg/bn/cond/NoOpNoOp^agg/bn/cond/switch_f

 agg/bn/cond/control_dependency_1Identityagg/bn/cond/switch_f^agg/bn/cond/NoOp*
T0
*'
_class
loc:@agg/bn/cond/switch_f*
_output_shapes
: 

agg/bn/cond/MergeMerge agg/bn/cond/control_dependency_1agg/bn/cond/control_dependency*
T0
*
N*
_output_shapes
: : 
_
agg/bn/cond_1/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
[
agg/bn/cond_1/switch_tIdentityagg/bn/cond_1/Switch:1*
T0
*
_output_shapes
: 
Y
agg/bn/cond_1/switch_fIdentityagg/bn/cond_1/Switch*
T0
*
_output_shapes
: 
Q
agg/bn/cond_1/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
}
agg/bn/cond_1/IdentityIdentityagg/bn/cond_1/Identity/Switch:1^agg/bn/cond/Merge*
T0*
_output_shapes	
:
®
agg/bn/cond_1/Identity/SwitchSwitchagg/bn/moments/Squeezeagg/bn/cond_1/pred_id*
T0*)
_class
loc:@agg/bn/moments/Squeeze*"
_output_shapes
::

agg/bn/cond_1/Identity_1Identity!agg/bn/cond_1/Identity_1/Switch:1^agg/bn/cond/Merge*
T0*
_output_shapes	
:
´
agg/bn/cond_1/Identity_1/SwitchSwitchagg/bn/moments/Squeeze_1agg/bn/cond_1/pred_id*
T0*+
_class!
loc:@agg/bn/moments/Squeeze_1*"
_output_shapes
::
ì
agg/bn/cond_1/Switch_1Switch;agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/readagg/bn/cond_1/pred_id*
T0*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*"
_output_shapes
::
ð
agg/bn/cond_1/Switch_2Switch=agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/readagg/bn/cond_1/pred_id*
T0*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*"
_output_shapes
::
}
agg/bn/cond_1/MergeMergeagg/bn/cond_1/Switch_1agg/bn/cond_1/Identity*
T0*
N*
_output_shapes
	:: 

agg/bn/cond_1/Merge_1Mergeagg/bn/cond_1/Switch_2agg/bn/cond_1/Identity_1*
T0*
N*
_output_shapes
	:: 
[
agg/bn/batchnorm/add/yConst*
valueB
 *o:*
dtype0*
_output_shapes
: 
p
agg/bn/batchnorm/addAddagg/bn/cond_1/Merge_1agg/bn/batchnorm/add/y*
T0*
_output_shapes	
:
[
agg/bn/batchnorm/RsqrtRsqrtagg/bn/batchnorm/add*
T0*
_output_shapes	
:
l
agg/bn/batchnorm/mulMulagg/bn/batchnorm/Rsqrtagg/bn/gamma/read*
T0*
_output_shapes	
:
r
agg/bn/batchnorm/mul_1Mulagg/BiasAddagg/bn/batchnorm/mul*
T0*'
_output_shapes
:2
n
agg/bn/batchnorm/mul_2Mulagg/bn/cond_1/Mergeagg/bn/batchnorm/mul*
T0*
_output_shapes	
:
k
agg/bn/batchnorm/subSubagg/bn/beta/readagg/bn/batchnorm/mul_2*
T0*
_output_shapes	
:
}
agg/bn/batchnorm/add_1Addagg/bn/batchnorm/mul_1agg/bn/batchnorm/sub*
T0*'
_output_shapes
:2
Z
agg/ReluReluagg/bn/batchnorm/add_1*
T0*'
_output_shapes
:2
©
avgpool/avgpoolAvgPoolagg/Relu*
ksize
2*
strides
*
paddingVALID*
data_formatNHWC*
T0*'
_output_shapes
:
i
Tile_5/multiplesConst*%
valueB"   2         *
dtype0*
_output_shapes
:
u
Tile_5Tileavgpool/avgpoolTile_5/multiples*
T0*

Tmultiples0*'
_output_shapes
:2
P
concat_10/axisConst*
value	B :*
dtype0*
_output_shapes
: 
~
	concat_10ConcatV2Tile_5agg/Reluconcat_10/axis*
N*
T0*

Tidx0*'
_output_shapes
:2
±
2seg/conv2/weights/Initializer/random_uniform/shapeConst*%
valueB"            *
dtype0*$
_class
loc:@seg/conv2/weights*
_output_shapes
:

0seg/conv2/weights/Initializer/random_uniform/minConst*
valueB
 *óµ½*
dtype0*$
_class
loc:@seg/conv2/weights*
_output_shapes
: 

0seg/conv2/weights/Initializer/random_uniform/maxConst*
valueB
 *óµ=*
dtype0*$
_class
loc:@seg/conv2/weights*
_output_shapes
: 
þ
:seg/conv2/weights/Initializer/random_uniform/RandomUniformRandomUniform2seg/conv2/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*$
_class
loc:@seg/conv2/weights*(
_output_shapes
:
â
0seg/conv2/weights/Initializer/random_uniform/subSub0seg/conv2/weights/Initializer/random_uniform/max0seg/conv2/weights/Initializer/random_uniform/min*
T0*$
_class
loc:@seg/conv2/weights*
_output_shapes
: 
þ
0seg/conv2/weights/Initializer/random_uniform/mulMul:seg/conv2/weights/Initializer/random_uniform/RandomUniform0seg/conv2/weights/Initializer/random_uniform/sub*
T0*$
_class
loc:@seg/conv2/weights*(
_output_shapes
:
ð
,seg/conv2/weights/Initializer/random_uniformAdd0seg/conv2/weights/Initializer/random_uniform/mul0seg/conv2/weights/Initializer/random_uniform/min*
T0*$
_class
loc:@seg/conv2/weights*(
_output_shapes
:
Î
seg/conv2/weights
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *$
_class
loc:@seg/conv2/weights*(
_output_shapes
:
ô
seg/conv2/weights/AssignAssignseg/conv2/weights,seg/conv2/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv2/weights*(
_output_shapes
:

seg/conv2/weights/readIdentityseg/conv2/weights"/device:CPU:0*
T0*$
_class
loc:@seg/conv2/weights*(
_output_shapes
:
Þ
seg/conv2/Conv2DConv2D	concat_10seg/conv2/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*'
_output_shapes
:2

"seg/conv2/biases/Initializer/ConstConst*
valueB*    *
dtype0*#
_class
loc:@seg/conv2/biases*
_output_shapes	
:
²
seg/conv2/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *#
_class
loc:@seg/conv2/biases*
_output_shapes	
:
Ú
seg/conv2/biases/AssignAssignseg/conv2/biases"seg/conv2/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv2/biases*
_output_shapes	
:

seg/conv2/biases/readIdentityseg/conv2/biases"/device:CPU:0*
T0*#
_class
loc:@seg/conv2/biases*
_output_shapes	
:

seg/conv2/BiasAddBiasAddseg/conv2/Conv2Dseg/conv2/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:2
¤
3seg/conv2/bn/beta/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*$
_class
loc:@seg/conv2/bn/beta*
_output_shapes
:

)seg/conv2/bn/beta/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*$
_class
loc:@seg/conv2/bn/beta*
_output_shapes
: 
é
#seg/conv2/bn/beta/Initializer/zerosFill3seg/conv2/bn/beta/Initializer/zeros/shape_as_tensor)seg/conv2/bn/beta/Initializer/zeros/Const*
T0*

index_type0*$
_class
loc:@seg/conv2/bn/beta*
_output_shapes	
:
´
seg/conv2/bn/beta
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *$
_class
loc:@seg/conv2/bn/beta*
_output_shapes	
:
Þ
seg/conv2/bn/beta/AssignAssignseg/conv2/bn/beta#seg/conv2/bn/beta/Initializer/zeros"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv2/bn/beta*
_output_shapes	
:

seg/conv2/bn/beta/readIdentityseg/conv2/bn/beta"/device:CPU:0*
T0*$
_class
loc:@seg/conv2/bn/beta*
_output_shapes	
:
¥
3seg/conv2/bn/gamma/Initializer/ones/shape_as_tensorConst*
valueB:*
dtype0*%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes
:

)seg/conv2/bn/gamma/Initializer/ones/ConstConst*
valueB
 *  ?*
dtype0*%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes
: 
ê
#seg/conv2/bn/gamma/Initializer/onesFill3seg/conv2/bn/gamma/Initializer/ones/shape_as_tensor)seg/conv2/bn/gamma/Initializer/ones/Const*
T0*

index_type0*%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes	
:
¶
seg/conv2/bn/gamma
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes	
:
á
seg/conv2/bn/gamma/AssignAssignseg/conv2/bn/gamma#seg/conv2/bn/gamma/Initializer/ones"/device:CPU:0*
T0*
validate_shape(*
use_locking(*%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes	
:

seg/conv2/bn/gamma/readIdentityseg/conv2/bn/gamma"/device:CPU:0*
T0*%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes	
:
¬
7seg/conv2/bn/pop_mean/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes
:

-seg/conv2/bn/pop_mean/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes
: 
ù
'seg/conv2/bn/pop_mean/Initializer/zerosFill7seg/conv2/bn/pop_mean/Initializer/zeros/shape_as_tensor-seg/conv2/bn/pop_mean/Initializer/zeros/Const*
T0*

index_type0*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes	
:
¼
seg/conv2/bn/pop_mean
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes	
:
î
seg/conv2/bn/pop_mean/AssignAssignseg/conv2/bn/pop_mean'seg/conv2/bn/pop_mean/Initializer/zeros"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes	
:

seg/conv2/bn/pop_mean/readIdentityseg/conv2/bn/pop_mean"/device:CPU:0*
T0*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes	
:
©
5seg/conv2/bn/pop_var/Initializer/ones/shape_as_tensorConst*
valueB:*
dtype0*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes
:

+seg/conv2/bn/pop_var/Initializer/ones/ConstConst*
valueB
 *  ?*
dtype0*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes
: 
ò
%seg/conv2/bn/pop_var/Initializer/onesFill5seg/conv2/bn/pop_var/Initializer/ones/shape_as_tensor+seg/conv2/bn/pop_var/Initializer/ones/Const*
T0*

index_type0*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes	
:
º
seg/conv2/bn/pop_var
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes	
:
é
seg/conv2/bn/pop_var/AssignAssignseg/conv2/bn/pop_var%seg/conv2/bn/pop_var/Initializer/ones"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes	
:

seg/conv2/bn/pop_var/readIdentityseg/conv2/bn/pop_var"/device:CPU:0*
T0*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes	
:
c
seg/conv2/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
c
seg/conv2/bn/cond/switch_tIdentityseg/conv2/bn/cond/Switch:1*
T0
*
_output_shapes
: 
a
seg/conv2/bn/cond/switch_fIdentityseg/conv2/bn/cond/Switch*
T0
*
_output_shapes
: 
U
seg/conv2/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
¢
0seg/conv2/bn/cond/moments/mean/reduction_indicesConst^seg/conv2/bn/cond/switch_t*!
valueB"          *
dtype0*
_output_shapes
:
Ð
seg/conv2/bn/cond/moments/meanMean'seg/conv2/bn/cond/moments/mean/Switch:10seg/conv2/bn/cond/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:
È
%seg/conv2/bn/cond/moments/mean/SwitchSwitchseg/conv2/BiasAddseg/conv2/bn/cond/pred_id*
T0*$
_class
loc:@seg/conv2/BiasAdd*:
_output_shapes(
&:2:2

&seg/conv2/bn/cond/moments/StopGradientStopGradientseg/conv2/bn/cond/moments/mean*
T0*'
_output_shapes
:
Ã
+seg/conv2/bn/cond/moments/SquaredDifferenceSquaredDifference'seg/conv2/bn/cond/moments/mean/Switch:1&seg/conv2/bn/cond/moments/StopGradient*
T0*'
_output_shapes
:2
¦
4seg/conv2/bn/cond/moments/variance/reduction_indicesConst^seg/conv2/bn/cond/switch_t*!
valueB"          *
dtype0*
_output_shapes
:
Ü
"seg/conv2/bn/cond/moments/varianceMean+seg/conv2/bn/cond/moments/SquaredDifference4seg/conv2/bn/cond/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:

!seg/conv2/bn/cond/moments/SqueezeSqueezeseg/conv2/bn/cond/moments/mean*
T0*
squeeze_dims
 *
_output_shapes	
:

#seg/conv2/bn/cond/moments/Squeeze_1Squeeze"seg/conv2/bn/cond/moments/variance*
T0*
squeeze_dims
 *
_output_shapes	
:
y
seg/conv2/bn/cond/mul/yConst^seg/conv2/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
{
seg/conv2/bn/cond/mulMulseg/conv2/bn/cond/mul/Switch:1seg/conv2/bn/cond/mul/y*
T0*
_output_shapes	
:
Ã
seg/conv2/bn/cond/mul/SwitchSwitchseg/conv2/bn/pop_mean/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*(
_class
loc:@seg/conv2/bn/pop_mean*"
_output_shapes
::
{
seg/conv2/bn/cond/mul_1/yConst^seg/conv2/bn/cond/switch_t*
valueB
 *ÍÌÌ=*
dtype0*
_output_shapes
: 

seg/conv2/bn/cond/mul_1Mul!seg/conv2/bn/cond/moments/Squeezeseg/conv2/bn/cond/mul_1/y*
T0*
_output_shapes	
:
r
seg/conv2/bn/cond/addAddseg/conv2/bn/cond/mulseg/conv2/bn/cond/mul_1*
T0*
_output_shapes	
:
ä
seg/conv2/bn/cond/AssignAssign!seg/conv2/bn/cond/Assign/Switch:1seg/conv2/bn/cond/add"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes	
:
Ä
seg/conv2/bn/cond/Assign/Switch	RefSwitchseg/conv2/bn/pop_meanseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*(
_class
loc:@seg/conv2/bn/pop_mean*"
_output_shapes
::
{
seg/conv2/bn/cond/mul_2/yConst^seg/conv2/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 

seg/conv2/bn/cond/mul_2Mul seg/conv2/bn/cond/mul_2/Switch:1seg/conv2/bn/cond/mul_2/y*
T0*
_output_shapes	
:
Ã
seg/conv2/bn/cond/mul_2/SwitchSwitchseg/conv2/bn/pop_var/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*'
_class
loc:@seg/conv2/bn/pop_var*"
_output_shapes
::
{
seg/conv2/bn/cond/mul_3/yConst^seg/conv2/bn/cond/switch_t*
valueB
 *ÍÌÌ=*
dtype0*
_output_shapes
: 

seg/conv2/bn/cond/mul_3Mul#seg/conv2/bn/cond/moments/Squeeze_1seg/conv2/bn/cond/mul_3/y*
T0*
_output_shapes	
:
v
seg/conv2/bn/cond/add_1Addseg/conv2/bn/cond/mul_2seg/conv2/bn/cond/mul_3*
T0*
_output_shapes	
:
é
seg/conv2/bn/cond/Assign_1Assign#seg/conv2/bn/cond/Assign_1/Switch:1seg/conv2/bn/cond/add_1"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes	
:
Ä
!seg/conv2/bn/cond/Assign_1/Switch	RefSwitchseg/conv2/bn/pop_varseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*'
_class
loc:@seg/conv2/bn/pop_var*"
_output_shapes
::
»
!seg/conv2/bn/cond/batchnorm/add/yConst^seg/conv2/bn/cond/Assign^seg/conv2/bn/cond/Assign_1^seg/conv2/bn/cond/switch_t*
valueB
 *o:*
dtype0*
_output_shapes
: 

seg/conv2/bn/cond/batchnorm/addAdd#seg/conv2/bn/cond/moments/Squeeze_1!seg/conv2/bn/cond/batchnorm/add/y*
T0*
_output_shapes	
:
q
!seg/conv2/bn/cond/batchnorm/RsqrtRsqrtseg/conv2/bn/cond/batchnorm/add*
T0*
_output_shapes	
:

seg/conv2/bn/cond/batchnorm/mulMul!seg/conv2/bn/cond/batchnorm/Rsqrt(seg/conv2/bn/cond/batchnorm/mul/Switch:1*
T0*
_output_shapes	
:
Ç
&seg/conv2/bn/cond/batchnorm/mul/SwitchSwitchseg/conv2/bn/gamma/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*%
_class
loc:@seg/conv2/bn/gamma*"
_output_shapes
::
¤
!seg/conv2/bn/cond/batchnorm/mul_1Mul'seg/conv2/bn/cond/moments/mean/Switch:1seg/conv2/bn/cond/batchnorm/mul*
T0*'
_output_shapes
:2

!seg/conv2/bn/cond/batchnorm/mul_2Mul!seg/conv2/bn/cond/moments/Squeezeseg/conv2/bn/cond/batchnorm/mul*
T0*
_output_shapes	
:

seg/conv2/bn/cond/batchnorm/subSub(seg/conv2/bn/cond/batchnorm/sub/Switch:1!seg/conv2/bn/cond/batchnorm/mul_2*
T0*
_output_shapes	
:
Å
&seg/conv2/bn/cond/batchnorm/sub/SwitchSwitchseg/conv2/bn/beta/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*$
_class
loc:@seg/conv2/bn/beta*"
_output_shapes
::

!seg/conv2/bn/cond/batchnorm/add_1Add!seg/conv2/bn/cond/batchnorm/mul_1seg/conv2/bn/cond/batchnorm/sub*
T0*'
_output_shapes
:2

#seg/conv2/bn/cond/batchnorm_1/add/yConst^seg/conv2/bn/cond/switch_f*
valueB
 *o:*
dtype0*
_output_shapes
: 

!seg/conv2/bn/cond/batchnorm_1/addAdd(seg/conv2/bn/cond/batchnorm_1/add/Switch#seg/conv2/bn/cond/batchnorm_1/add/y*
T0*
_output_shapes	
:
Í
(seg/conv2/bn/cond/batchnorm_1/add/SwitchSwitchseg/conv2/bn/pop_var/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*'
_class
loc:@seg/conv2/bn/pop_var*"
_output_shapes
::
u
#seg/conv2/bn/cond/batchnorm_1/RsqrtRsqrt!seg/conv2/bn/cond/batchnorm_1/add*
T0*
_output_shapes	
:

!seg/conv2/bn/cond/batchnorm_1/mulMul#seg/conv2/bn/cond/batchnorm_1/Rsqrt(seg/conv2/bn/cond/batchnorm_1/mul/Switch*
T0*
_output_shapes	
:
É
(seg/conv2/bn/cond/batchnorm_1/mul/SwitchSwitchseg/conv2/bn/gamma/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*%
_class
loc:@seg/conv2/bn/gamma*"
_output_shapes
::
«
#seg/conv2/bn/cond/batchnorm_1/mul_1Mul*seg/conv2/bn/cond/batchnorm_1/mul_1/Switch!seg/conv2/bn/cond/batchnorm_1/mul*
T0*'
_output_shapes
:2
Í
*seg/conv2/bn/cond/batchnorm_1/mul_1/SwitchSwitchseg/conv2/BiasAddseg/conv2/bn/cond/pred_id*
T0*$
_class
loc:@seg/conv2/BiasAdd*:
_output_shapes(
&:2:2

#seg/conv2/bn/cond/batchnorm_1/mul_2Mul*seg/conv2/bn/cond/batchnorm_1/mul_2/Switch!seg/conv2/bn/cond/batchnorm_1/mul*
T0*
_output_shapes	
:
Ñ
*seg/conv2/bn/cond/batchnorm_1/mul_2/SwitchSwitchseg/conv2/bn/pop_mean/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*(
_class
loc:@seg/conv2/bn/pop_mean*"
_output_shapes
::

!seg/conv2/bn/cond/batchnorm_1/subSub(seg/conv2/bn/cond/batchnorm_1/sub/Switch#seg/conv2/bn/cond/batchnorm_1/mul_2*
T0*
_output_shapes	
:
Ç
(seg/conv2/bn/cond/batchnorm_1/sub/SwitchSwitchseg/conv2/bn/beta/readseg/conv2/bn/cond/pred_id"/device:CPU:0*
T0*$
_class
loc:@seg/conv2/bn/beta*"
_output_shapes
::
¤
#seg/conv2/bn/cond/batchnorm_1/add_1Add#seg/conv2/bn/cond/batchnorm_1/mul_1!seg/conv2/bn/cond/batchnorm_1/sub*
T0*'
_output_shapes
:2
¥
seg/conv2/bn/cond/MergeMerge#seg/conv2/bn/cond/batchnorm_1/add_1!seg/conv2/bn/cond/batchnorm/add_1*
T0*
N*)
_output_shapes
:2: 
a
seg/conv2/ReluReluseg/conv2/bn/cond/Merge*
T0*'
_output_shapes
:2
^
seg/dp1/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
Y
seg/dp1/cond/switch_tIdentityseg/dp1/cond/Switch:1*
T0
*
_output_shapes
: 
W
seg/dp1/cond/switch_fIdentityseg/dp1/cond/Switch*
T0
*
_output_shapes
: 
P
seg/dp1/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
{
seg/dp1/cond/dropout/keep_probConst^seg/dp1/cond/switch_t*
valueB
 *?*
dtype0*
_output_shapes
: 

seg/dp1/cond/dropout/ShapeConst^seg/dp1/cond/switch_t*%
valueB"   2         *
dtype0*
_output_shapes
:

'seg/dp1/cond/dropout/random_uniform/minConst^seg/dp1/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

'seg/dp1/cond/dropout/random_uniform/maxConst^seg/dp1/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¶
1seg/dp1/cond/dropout/random_uniform/RandomUniformRandomUniformseg/dp1/cond/dropout/Shape*

seed *
seed2 *
dtype0*
T0*'
_output_shapes
:2
¡
'seg/dp1/cond/dropout/random_uniform/subSub'seg/dp1/cond/dropout/random_uniform/max'seg/dp1/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
¼
'seg/dp1/cond/dropout/random_uniform/mulMul1seg/dp1/cond/dropout/random_uniform/RandomUniform'seg/dp1/cond/dropout/random_uniform/sub*
T0*'
_output_shapes
:2
®
#seg/dp1/cond/dropout/random_uniformAdd'seg/dp1/cond/dropout/random_uniform/mul'seg/dp1/cond/dropout/random_uniform/min*
T0*'
_output_shapes
:2

seg/dp1/cond/dropout/addAddseg/dp1/cond/dropout/keep_prob#seg/dp1/cond/dropout/random_uniform*
T0*'
_output_shapes
:2
o
seg/dp1/cond/dropout/FloorFloorseg/dp1/cond/dropout/add*
T0*'
_output_shapes
:2

seg/dp1/cond/dropout/divRealDiv!seg/dp1/cond/dropout/div/Switch:1seg/dp1/cond/dropout/keep_prob*
T0*'
_output_shapes
:2
·
seg/dp1/cond/dropout/div/SwitchSwitchseg/conv2/Reluseg/dp1/cond/pred_id*
T0*!
_class
loc:@seg/conv2/Relu*:
_output_shapes(
&:2:2

seg/dp1/cond/dropout/mulMulseg/dp1/cond/dropout/divseg/dp1/cond/dropout/Floor*
T0*'
_output_shapes
:2
­
seg/dp1/cond/Switch_1Switchseg/conv2/Reluseg/dp1/cond/pred_id*
T0*!
_class
loc:@seg/conv2/Relu*:
_output_shapes(
&:2:2

seg/dp1/cond/MergeMergeseg/dp1/cond/Switch_1seg/dp1/cond/dropout/mul*
T0*
N*)
_output_shapes
:2: 
±
2seg/conv3/weights/Initializer/random_uniform/shapeConst*%
valueB"            *
dtype0*$
_class
loc:@seg/conv3/weights*
_output_shapes
:

0seg/conv3/weights/Initializer/random_uniform/minConst*
valueB
 *×³Ý½*
dtype0*$
_class
loc:@seg/conv3/weights*
_output_shapes
: 

0seg/conv3/weights/Initializer/random_uniform/maxConst*
valueB
 *×³Ý=*
dtype0*$
_class
loc:@seg/conv3/weights*
_output_shapes
: 
þ
:seg/conv3/weights/Initializer/random_uniform/RandomUniformRandomUniform2seg/conv3/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*$
_class
loc:@seg/conv3/weights*(
_output_shapes
:
â
0seg/conv3/weights/Initializer/random_uniform/subSub0seg/conv3/weights/Initializer/random_uniform/max0seg/conv3/weights/Initializer/random_uniform/min*
T0*$
_class
loc:@seg/conv3/weights*
_output_shapes
: 
þ
0seg/conv3/weights/Initializer/random_uniform/mulMul:seg/conv3/weights/Initializer/random_uniform/RandomUniform0seg/conv3/weights/Initializer/random_uniform/sub*
T0*$
_class
loc:@seg/conv3/weights*(
_output_shapes
:
ð
,seg/conv3/weights/Initializer/random_uniformAdd0seg/conv3/weights/Initializer/random_uniform/mul0seg/conv3/weights/Initializer/random_uniform/min*
T0*$
_class
loc:@seg/conv3/weights*(
_output_shapes
:
Î
seg/conv3/weights
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *$
_class
loc:@seg/conv3/weights*(
_output_shapes
:
ô
seg/conv3/weights/AssignAssignseg/conv3/weights,seg/conv3/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv3/weights*(
_output_shapes
:

seg/conv3/weights/readIdentityseg/conv3/weights"/device:CPU:0*
T0*$
_class
loc:@seg/conv3/weights*(
_output_shapes
:
ç
seg/conv3/Conv2DConv2Dseg/dp1/cond/Mergeseg/conv3/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*'
_output_shapes
:2

"seg/conv3/biases/Initializer/ConstConst*
valueB*    *
dtype0*#
_class
loc:@seg/conv3/biases*
_output_shapes	
:
²
seg/conv3/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *#
_class
loc:@seg/conv3/biases*
_output_shapes	
:
Ú
seg/conv3/biases/AssignAssignseg/conv3/biases"seg/conv3/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv3/biases*
_output_shapes	
:

seg/conv3/biases/readIdentityseg/conv3/biases"/device:CPU:0*
T0*#
_class
loc:@seg/conv3/biases*
_output_shapes	
:

seg/conv3/BiasAddBiasAddseg/conv3/Conv2Dseg/conv3/biases/read*
T0*
data_formatNHWC*'
_output_shapes
:2
¤
3seg/conv3/bn/beta/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*$
_class
loc:@seg/conv3/bn/beta*
_output_shapes
:

)seg/conv3/bn/beta/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*$
_class
loc:@seg/conv3/bn/beta*
_output_shapes
: 
é
#seg/conv3/bn/beta/Initializer/zerosFill3seg/conv3/bn/beta/Initializer/zeros/shape_as_tensor)seg/conv3/bn/beta/Initializer/zeros/Const*
T0*

index_type0*$
_class
loc:@seg/conv3/bn/beta*
_output_shapes	
:
´
seg/conv3/bn/beta
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *$
_class
loc:@seg/conv3/bn/beta*
_output_shapes	
:
Þ
seg/conv3/bn/beta/AssignAssignseg/conv3/bn/beta#seg/conv3/bn/beta/Initializer/zeros"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv3/bn/beta*
_output_shapes	
:

seg/conv3/bn/beta/readIdentityseg/conv3/bn/beta"/device:CPU:0*
T0*$
_class
loc:@seg/conv3/bn/beta*
_output_shapes	
:
¥
3seg/conv3/bn/gamma/Initializer/ones/shape_as_tensorConst*
valueB:*
dtype0*%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes
:

)seg/conv3/bn/gamma/Initializer/ones/ConstConst*
valueB
 *  ?*
dtype0*%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes
: 
ê
#seg/conv3/bn/gamma/Initializer/onesFill3seg/conv3/bn/gamma/Initializer/ones/shape_as_tensor)seg/conv3/bn/gamma/Initializer/ones/Const*
T0*

index_type0*%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes	
:
¶
seg/conv3/bn/gamma
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes	
:
á
seg/conv3/bn/gamma/AssignAssignseg/conv3/bn/gamma#seg/conv3/bn/gamma/Initializer/ones"/device:CPU:0*
T0*
validate_shape(*
use_locking(*%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes	
:

seg/conv3/bn/gamma/readIdentityseg/conv3/bn/gamma"/device:CPU:0*
T0*%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes	
:
¬
7seg/conv3/bn/pop_mean/Initializer/zeros/shape_as_tensorConst*
valueB:*
dtype0*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes
:

-seg/conv3/bn/pop_mean/Initializer/zeros/ConstConst*
valueB
 *    *
dtype0*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes
: 
ù
'seg/conv3/bn/pop_mean/Initializer/zerosFill7seg/conv3/bn/pop_mean/Initializer/zeros/shape_as_tensor-seg/conv3/bn/pop_mean/Initializer/zeros/Const*
T0*

index_type0*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes	
:
¼
seg/conv3/bn/pop_mean
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes	
:
î
seg/conv3/bn/pop_mean/AssignAssignseg/conv3/bn/pop_mean'seg/conv3/bn/pop_mean/Initializer/zeros"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes	
:

seg/conv3/bn/pop_mean/readIdentityseg/conv3/bn/pop_mean"/device:CPU:0*
T0*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes	
:
©
5seg/conv3/bn/pop_var/Initializer/ones/shape_as_tensorConst*
valueB:*
dtype0*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes
:

+seg/conv3/bn/pop_var/Initializer/ones/ConstConst*
valueB
 *  ?*
dtype0*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes
: 
ò
%seg/conv3/bn/pop_var/Initializer/onesFill5seg/conv3/bn/pop_var/Initializer/ones/shape_as_tensor+seg/conv3/bn/pop_var/Initializer/ones/Const*
T0*

index_type0*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes	
:
º
seg/conv3/bn/pop_var
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes	
:
é
seg/conv3/bn/pop_var/AssignAssignseg/conv3/bn/pop_var%seg/conv3/bn/pop_var/Initializer/ones"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes	
:

seg/conv3/bn/pop_var/readIdentityseg/conv3/bn/pop_var"/device:CPU:0*
T0*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes	
:
c
seg/conv3/bn/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
c
seg/conv3/bn/cond/switch_tIdentityseg/conv3/bn/cond/Switch:1*
T0
*
_output_shapes
: 
a
seg/conv3/bn/cond/switch_fIdentityseg/conv3/bn/cond/Switch*
T0
*
_output_shapes
: 
U
seg/conv3/bn/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
¢
0seg/conv3/bn/cond/moments/mean/reduction_indicesConst^seg/conv3/bn/cond/switch_t*!
valueB"          *
dtype0*
_output_shapes
:
Ð
seg/conv3/bn/cond/moments/meanMean'seg/conv3/bn/cond/moments/mean/Switch:10seg/conv3/bn/cond/moments/mean/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:
È
%seg/conv3/bn/cond/moments/mean/SwitchSwitchseg/conv3/BiasAddseg/conv3/bn/cond/pred_id*
T0*$
_class
loc:@seg/conv3/BiasAdd*:
_output_shapes(
&:2:2

&seg/conv3/bn/cond/moments/StopGradientStopGradientseg/conv3/bn/cond/moments/mean*
T0*'
_output_shapes
:
Ã
+seg/conv3/bn/cond/moments/SquaredDifferenceSquaredDifference'seg/conv3/bn/cond/moments/mean/Switch:1&seg/conv3/bn/cond/moments/StopGradient*
T0*'
_output_shapes
:2
¦
4seg/conv3/bn/cond/moments/variance/reduction_indicesConst^seg/conv3/bn/cond/switch_t*!
valueB"          *
dtype0*
_output_shapes
:
Ü
"seg/conv3/bn/cond/moments/varianceMean+seg/conv3/bn/cond/moments/SquaredDifference4seg/conv3/bn/cond/moments/variance/reduction_indices*
	keep_dims(*
T0*

Tidx0*'
_output_shapes
:

!seg/conv3/bn/cond/moments/SqueezeSqueezeseg/conv3/bn/cond/moments/mean*
T0*
squeeze_dims
 *
_output_shapes	
:

#seg/conv3/bn/cond/moments/Squeeze_1Squeeze"seg/conv3/bn/cond/moments/variance*
T0*
squeeze_dims
 *
_output_shapes	
:
y
seg/conv3/bn/cond/mul/yConst^seg/conv3/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 
{
seg/conv3/bn/cond/mulMulseg/conv3/bn/cond/mul/Switch:1seg/conv3/bn/cond/mul/y*
T0*
_output_shapes	
:
Ã
seg/conv3/bn/cond/mul/SwitchSwitchseg/conv3/bn/pop_mean/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*(
_class
loc:@seg/conv3/bn/pop_mean*"
_output_shapes
::
{
seg/conv3/bn/cond/mul_1/yConst^seg/conv3/bn/cond/switch_t*
valueB
 *ÍÌÌ=*
dtype0*
_output_shapes
: 

seg/conv3/bn/cond/mul_1Mul!seg/conv3/bn/cond/moments/Squeezeseg/conv3/bn/cond/mul_1/y*
T0*
_output_shapes	
:
r
seg/conv3/bn/cond/addAddseg/conv3/bn/cond/mulseg/conv3/bn/cond/mul_1*
T0*
_output_shapes	
:
ä
seg/conv3/bn/cond/AssignAssign!seg/conv3/bn/cond/Assign/Switch:1seg/conv3/bn/cond/add"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes	
:
Ä
seg/conv3/bn/cond/Assign/Switch	RefSwitchseg/conv3/bn/pop_meanseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*(
_class
loc:@seg/conv3/bn/pop_mean*"
_output_shapes
::
{
seg/conv3/bn/cond/mul_2/yConst^seg/conv3/bn/cond/switch_t*
valueB
 *fff?*
dtype0*
_output_shapes
: 

seg/conv3/bn/cond/mul_2Mul seg/conv3/bn/cond/mul_2/Switch:1seg/conv3/bn/cond/mul_2/y*
T0*
_output_shapes	
:
Ã
seg/conv3/bn/cond/mul_2/SwitchSwitchseg/conv3/bn/pop_var/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*'
_class
loc:@seg/conv3/bn/pop_var*"
_output_shapes
::
{
seg/conv3/bn/cond/mul_3/yConst^seg/conv3/bn/cond/switch_t*
valueB
 *ÍÌÌ=*
dtype0*
_output_shapes
: 

seg/conv3/bn/cond/mul_3Mul#seg/conv3/bn/cond/moments/Squeeze_1seg/conv3/bn/cond/mul_3/y*
T0*
_output_shapes	
:
v
seg/conv3/bn/cond/add_1Addseg/conv3/bn/cond/mul_2seg/conv3/bn/cond/mul_3*
T0*
_output_shapes	
:
é
seg/conv3/bn/cond/Assign_1Assign#seg/conv3/bn/cond/Assign_1/Switch:1seg/conv3/bn/cond/add_1"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes	
:
Ä
!seg/conv3/bn/cond/Assign_1/Switch	RefSwitchseg/conv3/bn/pop_varseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*'
_class
loc:@seg/conv3/bn/pop_var*"
_output_shapes
::
»
!seg/conv3/bn/cond/batchnorm/add/yConst^seg/conv3/bn/cond/Assign^seg/conv3/bn/cond/Assign_1^seg/conv3/bn/cond/switch_t*
valueB
 *o:*
dtype0*
_output_shapes
: 

seg/conv3/bn/cond/batchnorm/addAdd#seg/conv3/bn/cond/moments/Squeeze_1!seg/conv3/bn/cond/batchnorm/add/y*
T0*
_output_shapes	
:
q
!seg/conv3/bn/cond/batchnorm/RsqrtRsqrtseg/conv3/bn/cond/batchnorm/add*
T0*
_output_shapes	
:

seg/conv3/bn/cond/batchnorm/mulMul!seg/conv3/bn/cond/batchnorm/Rsqrt(seg/conv3/bn/cond/batchnorm/mul/Switch:1*
T0*
_output_shapes	
:
Ç
&seg/conv3/bn/cond/batchnorm/mul/SwitchSwitchseg/conv3/bn/gamma/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*%
_class
loc:@seg/conv3/bn/gamma*"
_output_shapes
::
¤
!seg/conv3/bn/cond/batchnorm/mul_1Mul'seg/conv3/bn/cond/moments/mean/Switch:1seg/conv3/bn/cond/batchnorm/mul*
T0*'
_output_shapes
:2

!seg/conv3/bn/cond/batchnorm/mul_2Mul!seg/conv3/bn/cond/moments/Squeezeseg/conv3/bn/cond/batchnorm/mul*
T0*
_output_shapes	
:

seg/conv3/bn/cond/batchnorm/subSub(seg/conv3/bn/cond/batchnorm/sub/Switch:1!seg/conv3/bn/cond/batchnorm/mul_2*
T0*
_output_shapes	
:
Å
&seg/conv3/bn/cond/batchnorm/sub/SwitchSwitchseg/conv3/bn/beta/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*$
_class
loc:@seg/conv3/bn/beta*"
_output_shapes
::

!seg/conv3/bn/cond/batchnorm/add_1Add!seg/conv3/bn/cond/batchnorm/mul_1seg/conv3/bn/cond/batchnorm/sub*
T0*'
_output_shapes
:2

#seg/conv3/bn/cond/batchnorm_1/add/yConst^seg/conv3/bn/cond/switch_f*
valueB
 *o:*
dtype0*
_output_shapes
: 

!seg/conv3/bn/cond/batchnorm_1/addAdd(seg/conv3/bn/cond/batchnorm_1/add/Switch#seg/conv3/bn/cond/batchnorm_1/add/y*
T0*
_output_shapes	
:
Í
(seg/conv3/bn/cond/batchnorm_1/add/SwitchSwitchseg/conv3/bn/pop_var/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*'
_class
loc:@seg/conv3/bn/pop_var*"
_output_shapes
::
u
#seg/conv3/bn/cond/batchnorm_1/RsqrtRsqrt!seg/conv3/bn/cond/batchnorm_1/add*
T0*
_output_shapes	
:

!seg/conv3/bn/cond/batchnorm_1/mulMul#seg/conv3/bn/cond/batchnorm_1/Rsqrt(seg/conv3/bn/cond/batchnorm_1/mul/Switch*
T0*
_output_shapes	
:
É
(seg/conv3/bn/cond/batchnorm_1/mul/SwitchSwitchseg/conv3/bn/gamma/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*%
_class
loc:@seg/conv3/bn/gamma*"
_output_shapes
::
«
#seg/conv3/bn/cond/batchnorm_1/mul_1Mul*seg/conv3/bn/cond/batchnorm_1/mul_1/Switch!seg/conv3/bn/cond/batchnorm_1/mul*
T0*'
_output_shapes
:2
Í
*seg/conv3/bn/cond/batchnorm_1/mul_1/SwitchSwitchseg/conv3/BiasAddseg/conv3/bn/cond/pred_id*
T0*$
_class
loc:@seg/conv3/BiasAdd*:
_output_shapes(
&:2:2

#seg/conv3/bn/cond/batchnorm_1/mul_2Mul*seg/conv3/bn/cond/batchnorm_1/mul_2/Switch!seg/conv3/bn/cond/batchnorm_1/mul*
T0*
_output_shapes	
:
Ñ
*seg/conv3/bn/cond/batchnorm_1/mul_2/SwitchSwitchseg/conv3/bn/pop_mean/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*(
_class
loc:@seg/conv3/bn/pop_mean*"
_output_shapes
::

!seg/conv3/bn/cond/batchnorm_1/subSub(seg/conv3/bn/cond/batchnorm_1/sub/Switch#seg/conv3/bn/cond/batchnorm_1/mul_2*
T0*
_output_shapes	
:
Ç
(seg/conv3/bn/cond/batchnorm_1/sub/SwitchSwitchseg/conv3/bn/beta/readseg/conv3/bn/cond/pred_id"/device:CPU:0*
T0*$
_class
loc:@seg/conv3/bn/beta*"
_output_shapes
::
¤
#seg/conv3/bn/cond/batchnorm_1/add_1Add#seg/conv3/bn/cond/batchnorm_1/mul_1!seg/conv3/bn/cond/batchnorm_1/sub*
T0*'
_output_shapes
:2
¥
seg/conv3/bn/cond/MergeMerge#seg/conv3/bn/cond/batchnorm_1/add_1!seg/conv3/bn/cond/batchnorm/add_1*
T0*
N*)
_output_shapes
:2: 
a
seg/conv3/ReluReluseg/conv3/bn/cond/Merge*
T0*'
_output_shapes
:2
^
seg/dp2/cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
Y
seg/dp2/cond/switch_tIdentityseg/dp2/cond/Switch:1*
T0
*
_output_shapes
: 
W
seg/dp2/cond/switch_fIdentityseg/dp2/cond/Switch*
T0
*
_output_shapes
: 
P
seg/dp2/cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
{
seg/dp2/cond/dropout/keep_probConst^seg/dp2/cond/switch_t*
valueB
 *?*
dtype0*
_output_shapes
: 

seg/dp2/cond/dropout/ShapeConst^seg/dp2/cond/switch_t*%
valueB"   2         *
dtype0*
_output_shapes
:

'seg/dp2/cond/dropout/random_uniform/minConst^seg/dp2/cond/switch_t*
valueB
 *    *
dtype0*
_output_shapes
: 

'seg/dp2/cond/dropout/random_uniform/maxConst^seg/dp2/cond/switch_t*
valueB
 *  ?*
dtype0*
_output_shapes
: 
¶
1seg/dp2/cond/dropout/random_uniform/RandomUniformRandomUniformseg/dp2/cond/dropout/Shape*

seed *
seed2 *
dtype0*
T0*'
_output_shapes
:2
¡
'seg/dp2/cond/dropout/random_uniform/subSub'seg/dp2/cond/dropout/random_uniform/max'seg/dp2/cond/dropout/random_uniform/min*
T0*
_output_shapes
: 
¼
'seg/dp2/cond/dropout/random_uniform/mulMul1seg/dp2/cond/dropout/random_uniform/RandomUniform'seg/dp2/cond/dropout/random_uniform/sub*
T0*'
_output_shapes
:2
®
#seg/dp2/cond/dropout/random_uniformAdd'seg/dp2/cond/dropout/random_uniform/mul'seg/dp2/cond/dropout/random_uniform/min*
T0*'
_output_shapes
:2

seg/dp2/cond/dropout/addAddseg/dp2/cond/dropout/keep_prob#seg/dp2/cond/dropout/random_uniform*
T0*'
_output_shapes
:2
o
seg/dp2/cond/dropout/FloorFloorseg/dp2/cond/dropout/add*
T0*'
_output_shapes
:2

seg/dp2/cond/dropout/divRealDiv!seg/dp2/cond/dropout/div/Switch:1seg/dp2/cond/dropout/keep_prob*
T0*'
_output_shapes
:2
·
seg/dp2/cond/dropout/div/SwitchSwitchseg/conv3/Reluseg/dp2/cond/pred_id*
T0*!
_class
loc:@seg/conv3/Relu*:
_output_shapes(
&:2:2

seg/dp2/cond/dropout/mulMulseg/dp2/cond/dropout/divseg/dp2/cond/dropout/Floor*
T0*'
_output_shapes
:2
­
seg/dp2/cond/Switch_1Switchseg/conv3/Reluseg/dp2/cond/pred_id*
T0*!
_class
loc:@seg/conv3/Relu*:
_output_shapes(
&:2:2

seg/dp2/cond/MergeMergeseg/dp2/cond/Switch_1seg/dp2/cond/dropout/mul*
T0*
N*)
_output_shapes
:2: 
±
2seg/conv5/weights/Initializer/random_uniform/shapeConst*%
valueB"            *
dtype0*$
_class
loc:@seg/conv5/weights*
_output_shapes
:

0seg/conv5/weights/Initializer/random_uniform/minConst*
valueB
 *(¾*
dtype0*$
_class
loc:@seg/conv5/weights*
_output_shapes
: 

0seg/conv5/weights/Initializer/random_uniform/maxConst*
valueB
 *(>*
dtype0*$
_class
loc:@seg/conv5/weights*
_output_shapes
: 
ý
:seg/conv5/weights/Initializer/random_uniform/RandomUniformRandomUniform2seg/conv5/weights/Initializer/random_uniform/shape*

seed *
seed2 *
dtype0*
T0*$
_class
loc:@seg/conv5/weights*'
_output_shapes
:
â
0seg/conv5/weights/Initializer/random_uniform/subSub0seg/conv5/weights/Initializer/random_uniform/max0seg/conv5/weights/Initializer/random_uniform/min*
T0*$
_class
loc:@seg/conv5/weights*
_output_shapes
: 
ý
0seg/conv5/weights/Initializer/random_uniform/mulMul:seg/conv5/weights/Initializer/random_uniform/RandomUniform0seg/conv5/weights/Initializer/random_uniform/sub*
T0*$
_class
loc:@seg/conv5/weights*'
_output_shapes
:
ï
,seg/conv5/weights/Initializer/random_uniformAdd0seg/conv5/weights/Initializer/random_uniform/mul0seg/conv5/weights/Initializer/random_uniform/min*
T0*$
_class
loc:@seg/conv5/weights*'
_output_shapes
:
Ì
seg/conv5/weights
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *$
_class
loc:@seg/conv5/weights*'
_output_shapes
:
ó
seg/conv5/weights/AssignAssignseg/conv5/weights,seg/conv5/weights/Initializer/random_uniform"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv5/weights*'
_output_shapes
:

seg/conv5/weights/readIdentityseg/conv5/weights"/device:CPU:0*
T0*$
_class
loc:@seg/conv5/weights*'
_output_shapes
:
æ
seg/conv5/Conv2DConv2Dseg/dp2/cond/Mergeseg/conv5/weights/read*
T0*
strides
*
use_cudnn_on_gpu(*
paddingVALID*
data_formatNHWC*
	dilations
*&
_output_shapes
:2

"seg/conv5/biases/Initializer/ConstConst*
valueB*    *
dtype0*#
_class
loc:@seg/conv5/biases*
_output_shapes
:
°
seg/conv5/biases
VariableV2"/device:CPU:0*
shape:*
dtype0*
	container *
shared_name *#
_class
loc:@seg/conv5/biases*
_output_shapes
:
Ù
seg/conv5/biases/AssignAssignseg/conv5/biases"seg/conv5/biases/Initializer/Const"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv5/biases*
_output_shapes
:

seg/conv5/biases/readIdentityseg/conv5/biases"/device:CPU:0*
T0*#
_class
loc:@seg/conv5/biases*
_output_shapes
:

seg/conv5/BiasAddBiasAddseg/conv5/Conv2Dseg/conv5/biases/read*
T0*
data_formatNHWC*&
_output_shapes
:2
V
cond/SwitchSwitchPlaceholder_3Placeholder_3*
T0
*
_output_shapes
: : 
I
cond/switch_tIdentitycond/Switch:1*
T0
*
_output_shapes
: 
G
cond/switch_fIdentitycond/Switch*
T0
*
_output_shapes
: 
H
cond/pred_idIdentityPlaceholder_3*
T0
*
_output_shapes
: 
¡
cond/Switch_1Switchseg/conv5/BiasAddcond/pred_id*
T0*$
_class
loc:@seg/conv5/BiasAdd*8
_output_shapes&
$:2:2
s

cond/ShapeConst^cond/switch_f*%
valueB"   2         *
dtype0*
_output_shapes
:
[
	cond/RankConst^cond/switch_f*
value	B :*
dtype0*
_output_shapes
: 
u
cond/Shape_1Const^cond/switch_f*%
valueB"   2         *
dtype0*
_output_shapes
:
\

cond/Sub/yConst^cond/switch_f*
value	B :*
dtype0*
_output_shapes
: 
G
cond/SubSub	cond/Rank
cond/Sub/y*
T0*
_output_shapes
: 
\
cond/Slice/beginPackcond/Sub*
N*
T0*

axis *
_output_shapes
:
i
cond/Slice/sizeConst^cond/switch_f*
valueB:*
dtype0*
_output_shapes
:
v

cond/SliceSlicecond/Shape_1cond/Slice/begincond/Slice/size*
T0*
Index0*
_output_shapes
:
w
cond/concat/values_0Const^cond/switch_f*
valueB:
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
:
b
cond/concat/axisConst^cond/switch_f*
value	B : *
dtype0*
_output_shapes
: 

cond/concatConcatV2cond/concat/values_0
cond/Slicecond/concat/axis*
N*
T0*

Tidx0*
_output_shapes
:
p
cond/ReshapeReshapecond/Reshape/Switchcond/concat*
T0*
Tshape0*
_output_shapes

:2
§
cond/Reshape/SwitchSwitchseg/conv5/BiasAddcond/pred_id*
T0*$
_class
loc:@seg/conv5/BiasAdd*8
_output_shapes&
$:2:2
N
cond/SoftmaxSoftmaxcond/Reshape*
T0*
_output_shapes

:2
r
cond/Reshape_1Reshapecond/Softmax
cond/Shape*
T0*
Tshape0*&
_output_shapes
:2
p

cond/MergeMergecond/Reshape_1cond/Switch_1:1*
T0*
N*(
_output_shapes
:2: 
e
Reshape_13/shapeConst*!
valueB"   2      *
dtype0*
_output_shapes
:
n

Reshape_13Reshape
cond/MergeReshape_13/shape*
T0*
Tshape0*"
_output_shapes
:2
z
)SparseSoftmaxCrossEntropyWithLogits/ShapeConst*
valueB"   2   *
dtype0*
_output_shapes
:

+SparseSoftmaxCrossEntropyWithLogits/Shape_1Const*!
valueB"   2      *
dtype0*
_output_shapes
:
j
(SparseSoftmaxCrossEntropyWithLogits/RankConst*
value	B :*
dtype0*
_output_shapes
: 
k
)SparseSoftmaxCrossEntropyWithLogits/sub/yConst*
value	B :*
dtype0*
_output_shapes
: 
¤
'SparseSoftmaxCrossEntropyWithLogits/subSub(SparseSoftmaxCrossEntropyWithLogits/Rank)SparseSoftmaxCrossEntropyWithLogits/sub/y*
T0*
_output_shapes
: 
k
)SparseSoftmaxCrossEntropyWithLogits/add/yConst*
value	B :*
dtype0*
_output_shapes
: 
£
'SparseSoftmaxCrossEntropyWithLogits/addAdd'SparseSoftmaxCrossEntropyWithLogits/sub)SparseSoftmaxCrossEntropyWithLogits/add/y*
T0*
_output_shapes
: 
¢
7SparseSoftmaxCrossEntropyWithLogits/strided_slice/stackPack'SparseSoftmaxCrossEntropyWithLogits/sub*
N*
T0*

axis *
_output_shapes
:
¤
9SparseSoftmaxCrossEntropyWithLogits/strided_slice/stack_1Pack'SparseSoftmaxCrossEntropyWithLogits/add*
N*
T0*

axis *
_output_shapes
:

9SparseSoftmaxCrossEntropyWithLogits/strided_slice/stack_2Const*
valueB:*
dtype0*
_output_shapes
:
¯
1SparseSoftmaxCrossEntropyWithLogits/strided_sliceStridedSlice+SparseSoftmaxCrossEntropyWithLogits/Shape_17SparseSoftmaxCrossEntropyWithLogits/strided_slice/stack9SparseSoftmaxCrossEntropyWithLogits/strided_slice/stack_19SparseSoftmaxCrossEntropyWithLogits/strided_slice/stack_2*
T0*
Index0*

begin_mask *
end_mask *
ellipsis_mask *
new_axis_mask *
shrink_axis_mask*
_output_shapes
: 
~
3SparseSoftmaxCrossEntropyWithLogits/Reshape/shape/0Const*
valueB :
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
: 
Û
1SparseSoftmaxCrossEntropyWithLogits/Reshape/shapePack3SparseSoftmaxCrossEntropyWithLogits/Reshape/shape/01SparseSoftmaxCrossEntropyWithLogits/strided_slice*
N*
T0*

axis *
_output_shapes
:
¾
+SparseSoftmaxCrossEntropyWithLogits/ReshapeReshape
Reshape_131SparseSoftmaxCrossEntropyWithLogits/Reshape/shape*
T0*
Tshape0*0
_output_shapes
:ÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿÿ

3SparseSoftmaxCrossEntropyWithLogits/Reshape_1/shapeConst*
valueB:
ÿÿÿÿÿÿÿÿÿ*
dtype0*
_output_shapes
:
¯
-SparseSoftmaxCrossEntropyWithLogits/Reshape_1ReshapePlaceholder_13SparseSoftmaxCrossEntropyWithLogits/Reshape_1/shape*
T0*
Tshape0*
_output_shapes
:2

GSparseSoftmaxCrossEntropyWithLogits/SparseSoftmaxCrossEntropyWithLogits#SparseSoftmaxCrossEntropyWithLogits+SparseSoftmaxCrossEntropyWithLogits/Reshape-SparseSoftmaxCrossEntropyWithLogits/Reshape_1*
T0*
Tlabels0*-
_output_shapes
:2:2ÿÿÿÿÿÿÿÿÿ
ã
-SparseSoftmaxCrossEntropyWithLogits/Reshape_2ReshapeGSparseSoftmaxCrossEntropyWithLogits/SparseSoftmaxCrossEntropyWithLogits)SparseSoftmaxCrossEntropyWithLogits/Shape*
T0*
Tshape0*
_output_shapes

:2
X
Mean/reduction_indicesConst*
value	B :*
dtype0*
_output_shapes
: 

MeanMean-SparseSoftmaxCrossEntropyWithLogits/Reshape_2Mean/reduction_indices*
	keep_dims( *
T0*

Tidx0*
_output_shapes
:
Q
Const_1Const*
valueB: *
dtype0*
_output_shapes
:
[
Mean_1MeanMeanConst_1*
	keep_dims( *
T0*

Tidx0*
_output_shapes
: 
R
ArgMax/dimensionConst*
value	B :*
dtype0*
_output_shapes
: 
v
ArgMaxArgMax
Reshape_13ArgMax/dimension*
T0*

Tidx0*
output_type0	*
_output_shapes

:2
P

save/ConstConst*
valueB Bmodel*
dtype0*
_output_shapes
: 
ã/
save/SaveV2/tensor_namesConst*/
value/B/{BBiasAdd/biasesBBiasAdd_1/biasesBBiasAdd_2/biasesBVariableB
agg/biasesB6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverageB8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageBagg/bn/betaBagg/bn/gammaBagg/weightsBgapnet00/biasesBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet00/weightsBgapnet01/biasesBgapnet01/bn/betaBgapnet01/bn/gammaB@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01/weightsBgapnet10/biasesBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/weightsBgapnet11/biasesBgapnet11/bn/betaBgapnet11/bn/gammaB@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11/weightsBglobal_expand/biasesBglobal_expand/bn/betaBglobal_expand/bn/gammaBJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverageBLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverageBglobal_expand/weightsBlayerfilter0_edgefea_0/biasesBlayerfilter0_edgefea_0/bn/betaBlayerfilter0_edgefea_0/bn/gammaB\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter0_edgefea_0/weightsB(layerfilter0_neib_att_conv_head_0/biasesB)layerfilter0_neib_att_conv_head_0/bn/betaB*layerfilter0_neib_att_conv_head_0/bn/gammaBrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_neib_att_conv_head_0/weightsB'layerfilter0_newfea_conv_head_0/bn/betaB(layerfilter0_newfea_conv_head_0/bn/gammaBnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter0_newfea_conv_head_0/weightsB(layerfilter0_self_att_conv_head_0/biasesB)layerfilter0_self_att_conv_head_0/bn/betaB*layerfilter0_self_att_conv_head_0/bn/gammaBrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_self_att_conv_head_0/weightsBlayerfilter1_edgefea_0/biasesBlayerfilter1_edgefea_0/bn/betaBlayerfilter1_edgefea_0/bn/gammaB\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_0/weightsBlayerfilter1_edgefea_1/biasesBlayerfilter1_edgefea_1/bn/betaBlayerfilter1_edgefea_1/bn/gammaB\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_1/weightsB(layerfilter1_neib_att_conv_head_0/biasesB)layerfilter1_neib_att_conv_head_0/bn/betaB*layerfilter1_neib_att_conv_head_0/bn/gammaBrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_0/weightsB(layerfilter1_neib_att_conv_head_1/biasesB)layerfilter1_neib_att_conv_head_1/bn/betaB*layerfilter1_neib_att_conv_head_1/bn/gammaBrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_1/weightsB'layerfilter1_newfea_conv_head_0/bn/betaB(layerfilter1_newfea_conv_head_0/bn/gammaBnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_0/weightsB'layerfilter1_newfea_conv_head_1/bn/betaB(layerfilter1_newfea_conv_head_1/bn/gammaBnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_1/weightsB(layerfilter1_self_att_conv_head_0/biasesB)layerfilter1_self_att_conv_head_0/bn/betaB*layerfilter1_self_att_conv_head_0/bn/gammaBrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_0/weightsB(layerfilter1_self_att_conv_head_1/biasesB)layerfilter1_self_att_conv_head_1/bn/betaB*layerfilter1_self_att_conv_head_1/bn/gammaBrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_1/weightsBseg/conv2/biasesBseg/conv2/bn/betaBseg/conv2/bn/gammaBseg/conv2/bn/pop_meanBseg/conv2/bn/pop_varBseg/conv2/weightsBseg/conv3/biasesBseg/conv3/bn/betaBseg/conv3/bn/gammaBseg/conv3/bn/pop_meanBseg/conv3/bn/pop_varBseg/conv3/weightsBseg/conv5/biasesBseg/conv5/weights*
dtype0*
_output_shapes
:{
Ü
save/SaveV2/shape_and_slicesConst*
valueBþ{B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:{
æ0
save/SaveV2SaveV2
save/Constsave/SaveV2/tensor_namessave/SaveV2/shape_and_slicesBiasAdd/biasesBiasAdd_1/biasesBiasAdd_2/biasesVariable
agg/biases6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageagg/bn/betaagg/bn/gammaagg/weightsgapnet00/biasesgapnet00/bn/betagapnet00/bn/gamma@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet00/weightsgapnet01/biasesgapnet01/bn/betagapnet01/bn/gamma@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverageBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet01/weightsgapnet10/biasesgapnet10/bn/betagapnet10/bn/gamma@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet10/weightsgapnet11/biasesgapnet11/bn/betagapnet11/bn/gamma@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverageBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet11/weightsglobal_expand/biasesglobal_expand/bn/betaglobal_expand/bn/gammaJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverageLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverageglobal_expand/weightslayerfilter0_edgefea_0/biaseslayerfilter0_edgefea_0/bn/betalayerfilter0_edgefea_0/bn/gamma\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0_edgefea_0/weights(layerfilter0_neib_att_conv_head_0/biases)layerfilter0_neib_att_conv_head_0/bn/beta*layerfilter0_neib_att_conv_head_0/bn/gammarlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter0_neib_att_conv_head_0/weights'layerfilter0_newfea_conv_head_0/bn/beta(layerfilter0_newfea_conv_head_0/bn/gammanlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage'layerfilter0_newfea_conv_head_0/weights(layerfilter0_self_att_conv_head_0/biases)layerfilter0_self_att_conv_head_0/bn/beta*layerfilter0_self_att_conv_head_0/bn/gammarlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter0_self_att_conv_head_0/weightslayerfilter1_edgefea_0/biaseslayerfilter1_edgefea_0/bn/betalayerfilter1_edgefea_0/bn/gamma\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_edgefea_0/weightslayerfilter1_edgefea_1/biaseslayerfilter1_edgefea_1/bn/betalayerfilter1_edgefea_1/bn/gamma\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_edgefea_1/weights(layerfilter1_neib_att_conv_head_0/biases)layerfilter1_neib_att_conv_head_0/bn/beta*layerfilter1_neib_att_conv_head_0/bn/gammarlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_neib_att_conv_head_0/weights(layerfilter1_neib_att_conv_head_1/biases)layerfilter1_neib_att_conv_head_1/bn/beta*layerfilter1_neib_att_conv_head_1/bn/gammarlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_neib_att_conv_head_1/weights'layerfilter1_newfea_conv_head_0/bn/beta(layerfilter1_newfea_conv_head_0/bn/gammanlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage'layerfilter1_newfea_conv_head_0/weights'layerfilter1_newfea_conv_head_1/bn/beta(layerfilter1_newfea_conv_head_1/bn/gammanlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage'layerfilter1_newfea_conv_head_1/weights(layerfilter1_self_att_conv_head_0/biases)layerfilter1_self_att_conv_head_0/bn/beta*layerfilter1_self_att_conv_head_0/bn/gammarlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_self_att_conv_head_0/weights(layerfilter1_self_att_conv_head_1/biases)layerfilter1_self_att_conv_head_1/bn/beta*layerfilter1_self_att_conv_head_1/bn/gammarlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_self_att_conv_head_1/weightsseg/conv2/biasesseg/conv2/bn/betaseg/conv2/bn/gammaseg/conv2/bn/pop_meanseg/conv2/bn/pop_varseg/conv2/weightsseg/conv3/biasesseg/conv3/bn/betaseg/conv3/bn/gammaseg/conv3/bn/pop_meanseg/conv3/bn/pop_varseg/conv3/weightsseg/conv5/biasesseg/conv5/weights*
dtypes
}2{
}
save/control_dependencyIdentity
save/Const^save/SaveV2*
T0*
_class
loc:@save/Const*
_output_shapes
: 
õ/
save/RestoreV2/tensor_namesConst"/device:CPU:0*/
value/B/{BBiasAdd/biasesBBiasAdd_1/biasesBBiasAdd_2/biasesBVariableB
agg/biasesB6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverageB8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageBagg/bn/betaBagg/bn/gammaBagg/weightsBgapnet00/biasesBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet00/weightsBgapnet01/biasesBgapnet01/bn/betaBgapnet01/bn/gammaB@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01/weightsBgapnet10/biasesBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/weightsBgapnet11/biasesBgapnet11/bn/betaBgapnet11/bn/gammaB@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11/weightsBglobal_expand/biasesBglobal_expand/bn/betaBglobal_expand/bn/gammaBJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverageBLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverageBglobal_expand/weightsBlayerfilter0_edgefea_0/biasesBlayerfilter0_edgefea_0/bn/betaBlayerfilter0_edgefea_0/bn/gammaB\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter0_edgefea_0/weightsB(layerfilter0_neib_att_conv_head_0/biasesB)layerfilter0_neib_att_conv_head_0/bn/betaB*layerfilter0_neib_att_conv_head_0/bn/gammaBrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_neib_att_conv_head_0/weightsB'layerfilter0_newfea_conv_head_0/bn/betaB(layerfilter0_newfea_conv_head_0/bn/gammaBnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter0_newfea_conv_head_0/weightsB(layerfilter0_self_att_conv_head_0/biasesB)layerfilter0_self_att_conv_head_0/bn/betaB*layerfilter0_self_att_conv_head_0/bn/gammaBrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_self_att_conv_head_0/weightsBlayerfilter1_edgefea_0/biasesBlayerfilter1_edgefea_0/bn/betaBlayerfilter1_edgefea_0/bn/gammaB\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_0/weightsBlayerfilter1_edgefea_1/biasesBlayerfilter1_edgefea_1/bn/betaBlayerfilter1_edgefea_1/bn/gammaB\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_1/weightsB(layerfilter1_neib_att_conv_head_0/biasesB)layerfilter1_neib_att_conv_head_0/bn/betaB*layerfilter1_neib_att_conv_head_0/bn/gammaBrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_0/weightsB(layerfilter1_neib_att_conv_head_1/biasesB)layerfilter1_neib_att_conv_head_1/bn/betaB*layerfilter1_neib_att_conv_head_1/bn/gammaBrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_1/weightsB'layerfilter1_newfea_conv_head_0/bn/betaB(layerfilter1_newfea_conv_head_0/bn/gammaBnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_0/weightsB'layerfilter1_newfea_conv_head_1/bn/betaB(layerfilter1_newfea_conv_head_1/bn/gammaBnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_1/weightsB(layerfilter1_self_att_conv_head_0/biasesB)layerfilter1_self_att_conv_head_0/bn/betaB*layerfilter1_self_att_conv_head_0/bn/gammaBrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_0/weightsB(layerfilter1_self_att_conv_head_1/biasesB)layerfilter1_self_att_conv_head_1/bn/betaB*layerfilter1_self_att_conv_head_1/bn/gammaBrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_1/weightsBseg/conv2/biasesBseg/conv2/bn/betaBseg/conv2/bn/gammaBseg/conv2/bn/pop_meanBseg/conv2/bn/pop_varBseg/conv2/weightsBseg/conv3/biasesBseg/conv3/bn/betaBseg/conv3/bn/gammaBseg/conv3/bn/pop_meanBseg/conv3/bn/pop_varBseg/conv3/weightsBseg/conv5/biasesBseg/conv5/weights*
dtype0*
_output_shapes
:{
î
save/RestoreV2/shape_and_slicesConst"/device:CPU:0*
valueBþ{B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:{

save/RestoreV2	RestoreV2
save/Constsave/RestoreV2/tensor_namessave/RestoreV2/shape_and_slices"/device:CPU:0*
dtypes
}2{*
_output_shapesï
ì:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
¦
save/AssignAssignBiasAdd/biasessave/RestoreV2*
T0*
validate_shape(*
use_locking(*!
_class
loc:@BiasAdd/biases*
_output_shapes
: 
®
save/Assign_1AssignBiasAdd_1/biasessave/RestoreV2:1*
T0*
validate_shape(*
use_locking(*#
_class
loc:@BiasAdd_1/biases*
_output_shapes
:@
®
save/Assign_2AssignBiasAdd_2/biasessave/RestoreV2:2*
T0*
validate_shape(*
use_locking(*#
_class
loc:@BiasAdd_2/biases*
_output_shapes
:@

save/Assign_3AssignVariablesave/RestoreV2:3*
T0*
validate_shape(*
use_locking(*
_class
loc:@Variable*
_output_shapes
: 
²
save/Assign_4Assign
agg/biasessave/RestoreV2:4"/device:CPU:0*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/biases*
_output_shapes	
:
û
save/Assign_5Assign6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:5*
T0*
validate_shape(*
use_locking(*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:
ÿ
save/Assign_6Assign8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:6*
T0*
validate_shape(*
use_locking(*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
¥
save/Assign_7Assignagg/bn/betasave/RestoreV2:7*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/bn/beta*
_output_shapes	
:
§
save/Assign_8Assignagg/bn/gammasave/RestoreV2:8*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/bn/gamma*
_output_shapes	
:
Á
save/Assign_9Assignagg/weightssave/RestoreV2:9"/device:CPU:0*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/weights*(
_output_shapes
:ð
½
save/Assign_10Assigngapnet00/biasessave/RestoreV2:10"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet00/biases*
_output_shapes
:@
°
save/Assign_11Assigngapnet00/bn/betasave/RestoreV2:11*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet00/bn/beta*
_output_shapes
:@
²
save/Assign_12Assigngapnet00/bn/gammasave/RestoreV2:12*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet00/bn/gamma*
_output_shapes
:@

save/Assign_13Assign@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:13*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

save/Assign_14AssignBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:14*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Ë
save/Assign_15Assigngapnet00/weightssave/RestoreV2:15"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@
¾
save/Assign_16Assigngapnet01/biasessave/RestoreV2:16"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet01/biases*
_output_shapes	
:
±
save/Assign_17Assigngapnet01/bn/betasave/RestoreV2:17*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet01/bn/beta*
_output_shapes	
:
³
save/Assign_18Assigngapnet01/bn/gammasave/RestoreV2:18*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet01/bn/gamma*
_output_shapes	
:

save/Assign_19Assign@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:19*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

save/Assign_20AssignBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:20*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
Ì
save/Assign_21Assigngapnet01/weightssave/RestoreV2:21"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet01/weights*'
_output_shapes
:@
¾
save/Assign_22Assigngapnet10/biasessave/RestoreV2:22"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet10/biases*
_output_shapes	
:
±
save/Assign_23Assigngapnet10/bn/betasave/RestoreV2:23*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet10/bn/beta*
_output_shapes	
:
³
save/Assign_24Assigngapnet10/bn/gammasave/RestoreV2:24*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet10/bn/gamma*
_output_shapes	
:

save/Assign_25Assign@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:25*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

save/Assign_26AssignBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:26*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
Í
save/Assign_27Assigngapnet10/weightssave/RestoreV2:27"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet10/weights*(
_output_shapes
:
¾
save/Assign_28Assigngapnet11/biasessave/RestoreV2:28"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet11/biases*
_output_shapes	
:
±
save/Assign_29Assigngapnet11/bn/betasave/RestoreV2:29*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet11/bn/beta*
_output_shapes	
:
³
save/Assign_30Assigngapnet11/bn/gammasave/RestoreV2:30*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet11/bn/gamma*
_output_shapes	
:

save/Assign_31Assign@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:31*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

save/Assign_32AssignBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:32*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
Í
save/Assign_33Assigngapnet11/weightssave/RestoreV2:33"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet11/weights*(
_output_shapes
:
Ç
save/Assign_34Assignglobal_expand/biasessave/RestoreV2:34"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@global_expand/biases*
_output_shapes
:
º
save/Assign_35Assignglobal_expand/bn/betasave/RestoreV2:35*
T0*
validate_shape(*
use_locking(*(
_class
loc:@global_expand/bn/beta*
_output_shapes
:
¼
save/Assign_36Assignglobal_expand/bn/gammasave/RestoreV2:36*
T0*
validate_shape(*
use_locking(*)
_class
loc:@global_expand/bn/gamma*
_output_shapes
:
¤
save/Assign_37AssignJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:37*
T0*
validate_shape(*
use_locking(*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¨
save/Assign_38AssignLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:38*
T0*
validate_shape(*
use_locking(*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Õ
save/Assign_39Assignglobal_expand/weightssave/RestoreV2:39"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@global_expand/weights*&
_output_shapes
:
Ù
save/Assign_40Assignlayerfilter0_edgefea_0/biasessave/RestoreV2:40"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter0_edgefea_0/biases*
_output_shapes
: 
Ì
save/Assign_41Assignlayerfilter0_edgefea_0/bn/betasave/RestoreV2:41*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter0_edgefea_0/bn/beta*
_output_shapes
: 
Î
save/Assign_42Assignlayerfilter0_edgefea_0/bn/gammasave/RestoreV2:42*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter0_edgefea_0/bn/gamma*
_output_shapes
: 
È
save/Assign_43Assign\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:43*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ì
save/Assign_44Assign^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:44*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ç
save/Assign_45Assignlayerfilter0_edgefea_0/weightssave/RestoreV2:45"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 
ï
save/Assign_46Assign(layerfilter0_neib_att_conv_head_0/biasessave/RestoreV2:46"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_neib_att_conv_head_0/biases*
_output_shapes
:
â
save/Assign_47Assign)layerfilter0_neib_att_conv_head_0/bn/betasave/RestoreV2:47*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/bn/beta*
_output_shapes
:
ä
save/Assign_48Assign*layerfilter0_neib_att_conv_head_0/bn/gammasave/RestoreV2:48*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter0_neib_att_conv_head_0/bn/gamma*
_output_shapes
:
õ
save/Assign_49Assignrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:49*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ù
save/Assign_50Assigntlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:50*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ý
save/Assign_51Assign)layerfilter0_neib_att_conv_head_0/weightssave/RestoreV2:51"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 
Þ
save/Assign_52Assign'layerfilter0_newfea_conv_head_0/bn/betasave/RestoreV2:52*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/bn/beta*
_output_shapes
: 
à
save/Assign_53Assign(layerfilter0_newfea_conv_head_0/bn/gammasave/RestoreV2:53*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_newfea_conv_head_0/bn/gamma*
_output_shapes
: 
í
save/Assign_54Assignnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:54*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
ñ
save/Assign_55Assignplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:55*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
ù
save/Assign_56Assign'layerfilter0_newfea_conv_head_0/weightssave/RestoreV2:56"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 
ï
save/Assign_57Assign(layerfilter0_self_att_conv_head_0/biasessave/RestoreV2:57"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_self_att_conv_head_0/biases*
_output_shapes
:
â
save/Assign_58Assign)layerfilter0_self_att_conv_head_0/bn/betasave/RestoreV2:58*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/bn/beta*
_output_shapes
:
ä
save/Assign_59Assign*layerfilter0_self_att_conv_head_0/bn/gammasave/RestoreV2:59*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter0_self_att_conv_head_0/bn/gamma*
_output_shapes
:
õ
save/Assign_60Assignrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:60*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ù
save/Assign_61Assigntlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:61*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ý
save/Assign_62Assign)layerfilter0_self_att_conv_head_0/weightssave/RestoreV2:62"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 
Ù
save/Assign_63Assignlayerfilter1_edgefea_0/biasessave/RestoreV2:63"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter1_edgefea_0/biases*
_output_shapes
:@
Ì
save/Assign_64Assignlayerfilter1_edgefea_0/bn/betasave/RestoreV2:64*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_0/bn/beta*
_output_shapes
:@
Î
save/Assign_65Assignlayerfilter1_edgefea_0/bn/gammasave/RestoreV2:65*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter1_edgefea_0/bn/gamma*
_output_shapes
:@
È
save/Assign_66Assign\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:66*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ì
save/Assign_67Assign^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:67*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
è
save/Assign_68Assignlayerfilter1_edgefea_0/weightssave/RestoreV2:68"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@
Ù
save/Assign_69Assignlayerfilter1_edgefea_1/biasessave/RestoreV2:69"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter1_edgefea_1/biases*
_output_shapes
:@
Ì
save/Assign_70Assignlayerfilter1_edgefea_1/bn/betasave/RestoreV2:70*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_1/bn/beta*
_output_shapes
:@
Î
save/Assign_71Assignlayerfilter1_edgefea_1/bn/gammasave/RestoreV2:71*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter1_edgefea_1/bn/gamma*
_output_shapes
:@
È
save/Assign_72Assign\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:72*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ì
save/Assign_73Assign^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:73*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
è
save/Assign_74Assignlayerfilter1_edgefea_1/weightssave/RestoreV2:74"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@
ï
save/Assign_75Assign(layerfilter1_neib_att_conv_head_0/biasessave/RestoreV2:75"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_0/biases*
_output_shapes
:
â
save/Assign_76Assign)layerfilter1_neib_att_conv_head_0/bn/betasave/RestoreV2:76*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/bn/beta*
_output_shapes
:
ä
save/Assign_77Assign*layerfilter1_neib_att_conv_head_0/bn/gammasave/RestoreV2:77*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_0/bn/gamma*
_output_shapes
:
õ
save/Assign_78Assignrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:78*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ù
save/Assign_79Assigntlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:79*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ý
save/Assign_80Assign)layerfilter1_neib_att_conv_head_0/weightssave/RestoreV2:80"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@
ï
save/Assign_81Assign(layerfilter1_neib_att_conv_head_1/biasessave/RestoreV2:81"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_1/biases*
_output_shapes
:
â
save/Assign_82Assign)layerfilter1_neib_att_conv_head_1/bn/betasave/RestoreV2:82*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/bn/beta*
_output_shapes
:
ä
save/Assign_83Assign*layerfilter1_neib_att_conv_head_1/bn/gammasave/RestoreV2:83*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_1/bn/gamma*
_output_shapes
:
õ
save/Assign_84Assignrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:84*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ù
save/Assign_85Assigntlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:85*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ý
save/Assign_86Assign)layerfilter1_neib_att_conv_head_1/weightssave/RestoreV2:86"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@
Þ
save/Assign_87Assign'layerfilter1_newfea_conv_head_0/bn/betasave/RestoreV2:87*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/bn/beta*
_output_shapes
:@
à
save/Assign_88Assign(layerfilter1_newfea_conv_head_0/bn/gammasave/RestoreV2:88*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_newfea_conv_head_0/bn/gamma*
_output_shapes
:@
í
save/Assign_89Assignnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:89*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ñ
save/Assign_90Assignplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:90*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
ú
save/Assign_91Assign'layerfilter1_newfea_conv_head_0/weightssave/RestoreV2:91"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@
Þ
save/Assign_92Assign'layerfilter1_newfea_conv_head_1/bn/betasave/RestoreV2:92*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/bn/beta*
_output_shapes
:@
à
save/Assign_93Assign(layerfilter1_newfea_conv_head_1/bn/gammasave/RestoreV2:93*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_newfea_conv_head_1/bn/gamma*
_output_shapes
:@
í
save/Assign_94Assignnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:94*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
ñ
save/Assign_95Assignplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:95*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
ú
save/Assign_96Assign'layerfilter1_newfea_conv_head_1/weightssave/RestoreV2:96"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@
ï
save/Assign_97Assign(layerfilter1_self_att_conv_head_0/biasessave/RestoreV2:97"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_self_att_conv_head_0/biases*
_output_shapes
:
â
save/Assign_98Assign)layerfilter1_self_att_conv_head_0/bn/betasave/RestoreV2:98*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/bn/beta*
_output_shapes
:
ä
save/Assign_99Assign*layerfilter1_self_att_conv_head_0/bn/gammasave/RestoreV2:99*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_self_att_conv_head_0/bn/gamma*
_output_shapes
:
÷
save/Assign_100Assignrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:100*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
û
save/Assign_101Assigntlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:101*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ÿ
save/Assign_102Assign)layerfilter1_self_att_conv_head_0/weightssave/RestoreV2:102"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@
ñ
save/Assign_103Assign(layerfilter1_self_att_conv_head_1/biasessave/RestoreV2:103"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_self_att_conv_head_1/biases*
_output_shapes
:
ä
save/Assign_104Assign)layerfilter1_self_att_conv_head_1/bn/betasave/RestoreV2:104*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/bn/beta*
_output_shapes
:
æ
save/Assign_105Assign*layerfilter1_self_att_conv_head_1/bn/gammasave/RestoreV2:105*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_self_att_conv_head_1/bn/gamma*
_output_shapes
:
÷
save/Assign_106Assignrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragesave/RestoreV2:106*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
û
save/Assign_107Assigntlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave/RestoreV2:107*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
ÿ
save/Assign_108Assign)layerfilter1_self_att_conv_head_1/weightssave/RestoreV2:108"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@
Â
save/Assign_109Assignseg/conv2/biasessave/RestoreV2:109"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv2/biases*
_output_shapes	
:
Ä
save/Assign_110Assignseg/conv2/bn/betasave/RestoreV2:110"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv2/bn/beta*
_output_shapes	
:
Æ
save/Assign_111Assignseg/conv2/bn/gammasave/RestoreV2:111"/device:CPU:0*
T0*
validate_shape(*
use_locking(*%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes	
:
Ì
save/Assign_112Assignseg/conv2/bn/pop_meansave/RestoreV2:112"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes	
:
Ê
save/Assign_113Assignseg/conv2/bn/pop_varsave/RestoreV2:113"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes	
:
Ñ
save/Assign_114Assignseg/conv2/weightssave/RestoreV2:114"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv2/weights*(
_output_shapes
:
Â
save/Assign_115Assignseg/conv3/biasessave/RestoreV2:115"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv3/biases*
_output_shapes	
:
Ä
save/Assign_116Assignseg/conv3/bn/betasave/RestoreV2:116"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv3/bn/beta*
_output_shapes	
:
Æ
save/Assign_117Assignseg/conv3/bn/gammasave/RestoreV2:117"/device:CPU:0*
T0*
validate_shape(*
use_locking(*%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes	
:
Ì
save/Assign_118Assignseg/conv3/bn/pop_meansave/RestoreV2:118"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes	
:
Ê
save/Assign_119Assignseg/conv3/bn/pop_varsave/RestoreV2:119"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes	
:
Ñ
save/Assign_120Assignseg/conv3/weightssave/RestoreV2:120"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv3/weights*(
_output_shapes
:
Á
save/Assign_121Assignseg/conv5/biasessave/RestoreV2:121"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv5/biases*
_output_shapes
:
Ð
save/Assign_122Assignseg/conv5/weightssave/RestoreV2:122"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv5/weights*'
_output_shapes
:
¥

save/restore_all/NoOpNoOp^save/Assign^save/Assign_1^save/Assign_2^save/Assign_3^save/Assign_5^save/Assign_6^save/Assign_7^save/Assign_8^save/Assign_11^save/Assign_12^save/Assign_13^save/Assign_14^save/Assign_17^save/Assign_18^save/Assign_19^save/Assign_20^save/Assign_23^save/Assign_24^save/Assign_25^save/Assign_26^save/Assign_29^save/Assign_30^save/Assign_31^save/Assign_32^save/Assign_35^save/Assign_36^save/Assign_37^save/Assign_38^save/Assign_41^save/Assign_42^save/Assign_43^save/Assign_44^save/Assign_47^save/Assign_48^save/Assign_49^save/Assign_50^save/Assign_52^save/Assign_53^save/Assign_54^save/Assign_55^save/Assign_58^save/Assign_59^save/Assign_60^save/Assign_61^save/Assign_64^save/Assign_65^save/Assign_66^save/Assign_67^save/Assign_70^save/Assign_71^save/Assign_72^save/Assign_73^save/Assign_76^save/Assign_77^save/Assign_78^save/Assign_79^save/Assign_82^save/Assign_83^save/Assign_84^save/Assign_85^save/Assign_87^save/Assign_88^save/Assign_89^save/Assign_90^save/Assign_92^save/Assign_93^save/Assign_94^save/Assign_95^save/Assign_98^save/Assign_99^save/Assign_100^save/Assign_101^save/Assign_104^save/Assign_105^save/Assign_106^save/Assign_107
Ü
save/restore_all/NoOp_1NoOp^save/Assign_4^save/Assign_9^save/Assign_10^save/Assign_15^save/Assign_16^save/Assign_21^save/Assign_22^save/Assign_27^save/Assign_28^save/Assign_33^save/Assign_34^save/Assign_39^save/Assign_40^save/Assign_45^save/Assign_46^save/Assign_51^save/Assign_56^save/Assign_57^save/Assign_62^save/Assign_63^save/Assign_68^save/Assign_69^save/Assign_74^save/Assign_75^save/Assign_80^save/Assign_81^save/Assign_86^save/Assign_91^save/Assign_96^save/Assign_97^save/Assign_102^save/Assign_103^save/Assign_108^save/Assign_109^save/Assign_110^save/Assign_111^save/Assign_112^save/Assign_113^save/Assign_114^save/Assign_115^save/Assign_116^save/Assign_117^save/Assign_118^save/Assign_119^save/Assign_120^save/Assign_121^save/Assign_122"/device:CPU:0
J
save/restore_allNoOp^save/restore_all/NoOp^save/restore_all/NoOp_1
R
save_1/ConstConst*
valueB Bmodel*
dtype0*
_output_shapes
: 

save_1/StringJoin/inputs_1Const*<
value3B1 B+_temp_36419754a984458b82984a9ea434af61/part*
dtype0*
_output_shapes
: 
{
save_1/StringJoin
StringJoinsave_1/Constsave_1/StringJoin/inputs_1*
N*
	separator *
_output_shapes
: 
S
save_1/num_shardsConst*
value	B :*
dtype0*
_output_shapes
: 
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
Õ%
save_1/SaveV2/tensor_namesConst"/device:CPU:0*÷$
valueí$Bê$LBBiasAdd/biasesBBiasAdd_1/biasesBBiasAdd_2/biasesBVariableB6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverageB8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageBagg/bn/betaBagg/bn/gammaBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01/bn/betaBgapnet01/bn/gammaB@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11/bn/betaBgapnet11/bn/gammaB@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverageBglobal_expand/bn/betaBglobal_expand/bn/gammaBJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverageBLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter0_edgefea_0/bn/betaBlayerfilter0_edgefea_0/bn/gammaB\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_neib_att_conv_head_0/bn/betaB*layerfilter0_neib_att_conv_head_0/bn/gammaBrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter0_newfea_conv_head_0/bn/betaB(layerfilter0_newfea_conv_head_0/bn/gammaBnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_self_att_conv_head_0/bn/betaB*layerfilter0_self_att_conv_head_0/bn/gammaBrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_0/bn/betaBlayerfilter1_edgefea_0/bn/gammaB\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_1/bn/betaBlayerfilter1_edgefea_1/bn/gammaB\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_0/bn/betaB*layerfilter1_neib_att_conv_head_0/bn/gammaBrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_1/bn/betaB*layerfilter1_neib_att_conv_head_1/bn/gammaBrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_0/bn/betaB(layerfilter1_newfea_conv_head_0/bn/gammaBnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_1/bn/betaB(layerfilter1_newfea_conv_head_1/bn/gammaBnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_0/bn/betaB*layerfilter1_self_att_conv_head_0/bn/gammaBrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_1/bn/betaB*layerfilter1_self_att_conv_head_1/bn/gammaBrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:L

save_1/SaveV2/shape_and_slicesConst"/device:CPU:0*­
value£B LB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:L
¸&
save_1/SaveV2SaveV2save_1/ShardedFilenamesave_1/SaveV2/tensor_namessave_1/SaveV2/shape_and_slicesBiasAdd/biasesBiasAdd_1/biasesBiasAdd_2/biasesVariable6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageagg/bn/betaagg/bn/gammagapnet00/bn/betagapnet00/bn/gamma@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet01/bn/betagapnet01/bn/gamma@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverageBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet10/bn/betagapnet10/bn/gamma@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragegapnet11/bn/betagapnet11/bn/gamma@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverageBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverageglobal_expand/bn/betaglobal_expand/bn/gammaJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverageLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter0_edgefea_0/bn/betalayerfilter0_edgefea_0/bn/gamma\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter0_neib_att_conv_head_0/bn/beta*layerfilter0_neib_att_conv_head_0/bn/gammarlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage'layerfilter0_newfea_conv_head_0/bn/beta(layerfilter0_newfea_conv_head_0/bn/gammanlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter0_self_att_conv_head_0/bn/beta*layerfilter0_self_att_conv_head_0/bn/gammarlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_edgefea_0/bn/betalayerfilter1_edgefea_0/bn/gamma\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragelayerfilter1_edgefea_1/bn/betalayerfilter1_edgefea_1/bn/gamma\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_neib_att_conv_head_0/bn/beta*layerfilter1_neib_att_conv_head_0/bn/gammarlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_neib_att_conv_head_1/bn/beta*layerfilter1_neib_att_conv_head_1/bn/gammarlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage'layerfilter1_newfea_conv_head_0/bn/beta(layerfilter1_newfea_conv_head_0/bn/gammanlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage'layerfilter1_newfea_conv_head_1/bn/beta(layerfilter1_newfea_conv_head_1/bn/gammanlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_self_att_conv_head_0/bn/beta*layerfilter1_self_att_conv_head_0/bn/gammarlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage)layerfilter1_self_att_conv_head_1/bn/beta*layerfilter1_self_att_conv_head_1/bn/gammarlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragetlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage"/device:CPU:0*Z
dtypesP
N2L
¨
save_1/control_dependencyIdentitysave_1/ShardedFilename^save_1/SaveV2"/device:CPU:0*
T0*)
_class
loc:@save_1/ShardedFilename*
_output_shapes
: 
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

save_1/SaveV2_1/tensor_namesConst"/device:CPU:0*´

valueª
B§
/B
agg/biasesBagg/weightsBgapnet00/biasesBgapnet00/weightsBgapnet01/biasesBgapnet01/weightsBgapnet10/biasesBgapnet10/weightsBgapnet11/biasesBgapnet11/weightsBglobal_expand/biasesBglobal_expand/weightsBlayerfilter0_edgefea_0/biasesBlayerfilter0_edgefea_0/weightsB(layerfilter0_neib_att_conv_head_0/biasesB)layerfilter0_neib_att_conv_head_0/weightsB'layerfilter0_newfea_conv_head_0/weightsB(layerfilter0_self_att_conv_head_0/biasesB)layerfilter0_self_att_conv_head_0/weightsBlayerfilter1_edgefea_0/biasesBlayerfilter1_edgefea_0/weightsBlayerfilter1_edgefea_1/biasesBlayerfilter1_edgefea_1/weightsB(layerfilter1_neib_att_conv_head_0/biasesB)layerfilter1_neib_att_conv_head_0/weightsB(layerfilter1_neib_att_conv_head_1/biasesB)layerfilter1_neib_att_conv_head_1/weightsB'layerfilter1_newfea_conv_head_0/weightsB'layerfilter1_newfea_conv_head_1/weightsB(layerfilter1_self_att_conv_head_0/biasesB)layerfilter1_self_att_conv_head_0/weightsB(layerfilter1_self_att_conv_head_1/biasesB)layerfilter1_self_att_conv_head_1/weightsBseg/conv2/biasesBseg/conv2/bn/betaBseg/conv2/bn/gammaBseg/conv2/bn/pop_meanBseg/conv2/bn/pop_varBseg/conv2/weightsBseg/conv3/biasesBseg/conv3/bn/betaBseg/conv3/bn/gammaBseg/conv3/bn/pop_meanBseg/conv3/bn/pop_varBseg/conv3/weightsBseg/conv5/biasesBseg/conv5/weights*
dtype0*
_output_shapes
:/
Ô
 save_1/SaveV2_1/shape_and_slicesConst"/device:CPU:0*q
valuehBf/B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:/
à
save_1/SaveV2_1SaveV2save_1/ShardedFilename_1save_1/SaveV2_1/tensor_names save_1/SaveV2_1/shape_and_slices
agg/biasesagg/weightsgapnet00/biasesgapnet00/weightsgapnet01/biasesgapnet01/weightsgapnet10/biasesgapnet10/weightsgapnet11/biasesgapnet11/weightsglobal_expand/biasesglobal_expand/weightslayerfilter0_edgefea_0/biaseslayerfilter0_edgefea_0/weights(layerfilter0_neib_att_conv_head_0/biases)layerfilter0_neib_att_conv_head_0/weights'layerfilter0_newfea_conv_head_0/weights(layerfilter0_self_att_conv_head_0/biases)layerfilter0_self_att_conv_head_0/weightslayerfilter1_edgefea_0/biaseslayerfilter1_edgefea_0/weightslayerfilter1_edgefea_1/biaseslayerfilter1_edgefea_1/weights(layerfilter1_neib_att_conv_head_0/biases)layerfilter1_neib_att_conv_head_0/weights(layerfilter1_neib_att_conv_head_1/biases)layerfilter1_neib_att_conv_head_1/weights'layerfilter1_newfea_conv_head_0/weights'layerfilter1_newfea_conv_head_1/weights(layerfilter1_self_att_conv_head_0/biases)layerfilter1_self_att_conv_head_0/weights(layerfilter1_self_att_conv_head_1/biases)layerfilter1_self_att_conv_head_1/weightsseg/conv2/biasesseg/conv2/bn/betaseg/conv2/bn/gammaseg/conv2/bn/pop_meanseg/conv2/bn/pop_varseg/conv2/weightsseg/conv3/biasesseg/conv3/bn/betaseg/conv3/bn/gammaseg/conv3/bn/pop_meanseg/conv3/bn/pop_varseg/conv3/weightsseg/conv5/biasesseg/conv5/weights"/device:CPU:0*=
dtypes3
12/
°
save_1/control_dependency_1Identitysave_1/ShardedFilename_1^save_1/SaveV2_1"/device:CPU:0*
T0*+
_class!
loc:@save_1/ShardedFilename_1*
_output_shapes
: 
ê
-save_1/MergeV2Checkpoints/checkpoint_prefixesPacksave_1/ShardedFilenamesave_1/ShardedFilename_1^save_1/control_dependency^save_1/control_dependency_1"/device:CPU:0*
N*
T0*

axis *
_output_shapes
:

save_1/MergeV2CheckpointsMergeV2Checkpoints-save_1/MergeV2Checkpoints/checkpoint_prefixessave_1/Const"/device:CPU:0*
delete_old_dirs(
¯
save_1/IdentityIdentitysave_1/Const^save_1/control_dependency^save_1/control_dependency_1^save_1/MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 
Ø%
save_1/RestoreV2/tensor_namesConst"/device:CPU:0*÷$
valueí$Bê$LBBiasAdd/biasesBBiasAdd_1/biasesBBiasAdd_2/biasesBVariableB6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverageB8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverageBagg/bn/betaBagg/bn/gammaBgapnet00/bn/betaBgapnet00/bn/gammaB@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet01/bn/betaBgapnet01/bn/gammaB@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet10/bn/betaBgapnet10/bn/gammaB@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverageBgapnet11/bn/betaBgapnet11/bn/gammaB@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverageBBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverageBglobal_expand/bn/betaBglobal_expand/bn/gammaBJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverageBLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter0_edgefea_0/bn/betaBlayerfilter0_edgefea_0/bn/gammaB\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_neib_att_conv_head_0/bn/betaB*layerfilter0_neib_att_conv_head_0/bn/gammaBrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter0_newfea_conv_head_0/bn/betaB(layerfilter0_newfea_conv_head_0/bn/gammaBnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter0_self_att_conv_head_0/bn/betaB*layerfilter0_self_att_conv_head_0/bn/gammaBrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_0/bn/betaBlayerfilter1_edgefea_0/bn/gammaB\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverageBlayerfilter1_edgefea_1/bn/betaBlayerfilter1_edgefea_1/bn/gammaB\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverageB^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_0/bn/betaB*layerfilter1_neib_att_conv_head_0/bn/gammaBrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_neib_att_conv_head_1/bn/betaB*layerfilter1_neib_att_conv_head_1/bn/gammaBrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_0/bn/betaB(layerfilter1_newfea_conv_head_0/bn/gammaBnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB'layerfilter1_newfea_conv_head_1/bn/betaB(layerfilter1_newfea_conv_head_1/bn/gammaBnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_0/bn/betaB*layerfilter1_self_att_conv_head_0/bn/gammaBrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverageB)layerfilter1_self_att_conv_head_1/bn/betaB*layerfilter1_self_att_conv_head_1/bn/gammaBrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverageBtlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
dtype0*
_output_shapes
:L

!save_1/RestoreV2/shape_and_slicesConst"/device:CPU:0*­
value£B LB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:L
¡
save_1/RestoreV2	RestoreV2save_1/Constsave_1/RestoreV2/tensor_names!save_1/RestoreV2/shape_and_slices"/device:CPU:0*Z
dtypesP
N2L*Æ
_output_shapes³
°::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ª
save_1/AssignAssignBiasAdd/biasessave_1/RestoreV2*
T0*
validate_shape(*
use_locking(*!
_class
loc:@BiasAdd/biases*
_output_shapes
: 
²
save_1/Assign_1AssignBiasAdd_1/biasessave_1/RestoreV2:1*
T0*
validate_shape(*
use_locking(*#
_class
loc:@BiasAdd_1/biases*
_output_shapes
:@
²
save_1/Assign_2AssignBiasAdd_2/biasessave_1/RestoreV2:2*
T0*
validate_shape(*
use_locking(*#
_class
loc:@BiasAdd_2/biases*
_output_shapes
:@

save_1/Assign_3AssignVariablesave_1/RestoreV2:3*
T0*
validate_shape(*
use_locking(*
_class
loc:@Variable*
_output_shapes
: 
ÿ
save_1/Assign_4Assign6agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:4*
T0*
validate_shape(*
use_locking(*I
_class?
=;loc:@agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

save_1/Assign_5Assign8agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:5*
T0*
validate_shape(*
use_locking(*K
_classA
?=loc:@agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
©
save_1/Assign_6Assignagg/bn/betasave_1/RestoreV2:6*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/bn/beta*
_output_shapes	
:
«
save_1/Assign_7Assignagg/bn/gammasave_1/RestoreV2:7*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/bn/gamma*
_output_shapes	
:
²
save_1/Assign_8Assigngapnet00/bn/betasave_1/RestoreV2:8*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet00/bn/beta*
_output_shapes
:@
´
save_1/Assign_9Assigngapnet00/bn/gammasave_1/RestoreV2:9*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet00/bn/gamma*
_output_shapes
:@

save_1/Assign_10Assign@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:10*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@

save_1/Assign_11AssignBgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:11*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
µ
save_1/Assign_12Assigngapnet01/bn/betasave_1/RestoreV2:12*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet01/bn/beta*
_output_shapes	
:
·
save_1/Assign_13Assigngapnet01/bn/gammasave_1/RestoreV2:13*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet01/bn/gamma*
_output_shapes	
:

save_1/Assign_14Assign@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:14*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

save_1/Assign_15AssignBgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:15*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
µ
save_1/Assign_16Assigngapnet10/bn/betasave_1/RestoreV2:16*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet10/bn/beta*
_output_shapes	
:
·
save_1/Assign_17Assigngapnet10/bn/gammasave_1/RestoreV2:17*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet10/bn/gamma*
_output_shapes	
:

save_1/Assign_18Assign@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:18*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

save_1/Assign_19AssignBgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:19*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
µ
save_1/Assign_20Assigngapnet11/bn/betasave_1/RestoreV2:20*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet11/bn/beta*
_output_shapes	
:
·
save_1/Assign_21Assigngapnet11/bn/gammasave_1/RestoreV2:21*
T0*
validate_shape(*
use_locking(*$
_class
loc:@gapnet11/bn/gamma*
_output_shapes	
:

save_1/Assign_22Assign@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:22*
T0*
validate_shape(*
use_locking(*S
_classI
GEloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes	
:

save_1/Assign_23AssignBgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:23*
T0*
validate_shape(*
use_locking(*U
_classK
IGloc:@gapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes	
:
¾
save_1/Assign_24Assignglobal_expand/bn/betasave_1/RestoreV2:24*
T0*
validate_shape(*
use_locking(*(
_class
loc:@global_expand/bn/beta*
_output_shapes
:
À
save_1/Assign_25Assignglobal_expand/bn/gammasave_1/RestoreV2:25*
T0*
validate_shape(*
use_locking(*)
_class
loc:@global_expand/bn/gamma*
_output_shapes
:
¨
save_1/Assign_26AssignJglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:26*
T0*
validate_shape(*
use_locking(*]
_classS
QOloc:@global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
¬
save_1/Assign_27AssignLglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:27*
T0*
validate_shape(*
use_locking(*_
_classU
SQloc:@global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ð
save_1/Assign_28Assignlayerfilter0_edgefea_0/bn/betasave_1/RestoreV2:28*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter0_edgefea_0/bn/beta*
_output_shapes
: 
Ò
save_1/Assign_29Assignlayerfilter0_edgefea_0/bn/gammasave_1/RestoreV2:29*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter0_edgefea_0/bn/gamma*
_output_shapes
: 
Ì
save_1/Assign_30Assign\layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:30*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
Ð
save_1/Assign_31Assign^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:31*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
æ
save_1/Assign_32Assign)layerfilter0_neib_att_conv_head_0/bn/betasave_1/RestoreV2:32*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/bn/beta*
_output_shapes
:
è
save_1/Assign_33Assign*layerfilter0_neib_att_conv_head_0/bn/gammasave_1/RestoreV2:33*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter0_neib_att_conv_head_0/bn/gamma*
_output_shapes
:
ù
save_1/Assign_34Assignrlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:34*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ý
save_1/Assign_35Assigntlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:35*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
â
save_1/Assign_36Assign'layerfilter0_newfea_conv_head_0/bn/betasave_1/RestoreV2:36*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/bn/beta*
_output_shapes
: 
ä
save_1/Assign_37Assign(layerfilter0_newfea_conv_head_0/bn/gammasave_1/RestoreV2:37*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_newfea_conv_head_0/bn/gamma*
_output_shapes
: 
ñ
save_1/Assign_38Assignnlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:38*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
: 
õ
save_1/Assign_39Assignplayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:39*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
: 
æ
save_1/Assign_40Assign)layerfilter0_self_att_conv_head_0/bn/betasave_1/RestoreV2:40*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/bn/beta*
_output_shapes
:
è
save_1/Assign_41Assign*layerfilter0_self_att_conv_head_0/bn/gammasave_1/RestoreV2:41*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter0_self_att_conv_head_0/bn/gamma*
_output_shapes
:
ù
save_1/Assign_42Assignrlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:42*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ý
save_1/Assign_43Assigntlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:43*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
Ð
save_1/Assign_44Assignlayerfilter1_edgefea_0/bn/betasave_1/RestoreV2:44*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_0/bn/beta*
_output_shapes
:@
Ò
save_1/Assign_45Assignlayerfilter1_edgefea_0/bn/gammasave_1/RestoreV2:45*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter1_edgefea_0/bn/gamma*
_output_shapes
:@
Ì
save_1/Assign_46Assign\layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:46*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ð
save_1/Assign_47Assign^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:47*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
Ð
save_1/Assign_48Assignlayerfilter1_edgefea_1/bn/betasave_1/RestoreV2:48*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_1/bn/beta*
_output_shapes
:@
Ò
save_1/Assign_49Assignlayerfilter1_edgefea_1/bn/gammasave_1/RestoreV2:49*
T0*
validate_shape(*
use_locking(*2
_class(
&$loc:@layerfilter1_edgefea_1/bn/gamma*
_output_shapes
:@
Ì
save_1/Assign_50Assign\layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:50*
T0*
validate_shape(*
use_locking(*o
_classe
caloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
Ð
save_1/Assign_51Assign^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:51*
T0*
validate_shape(*
use_locking(*q
_classg
ecloc:@layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
æ
save_1/Assign_52Assign)layerfilter1_neib_att_conv_head_0/bn/betasave_1/RestoreV2:52*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/bn/beta*
_output_shapes
:
è
save_1/Assign_53Assign*layerfilter1_neib_att_conv_head_0/bn/gammasave_1/RestoreV2:53*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_0/bn/gamma*
_output_shapes
:
ù
save_1/Assign_54Assignrlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:54*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ý
save_1/Assign_55Assigntlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:55*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
æ
save_1/Assign_56Assign)layerfilter1_neib_att_conv_head_1/bn/betasave_1/RestoreV2:56*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/bn/beta*
_output_shapes
:
è
save_1/Assign_57Assign*layerfilter1_neib_att_conv_head_1/bn/gammasave_1/RestoreV2:57*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_neib_att_conv_head_1/bn/gamma*
_output_shapes
:
ù
save_1/Assign_58Assignrlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:58*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ý
save_1/Assign_59Assigntlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:59*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
â
save_1/Assign_60Assign'layerfilter1_newfea_conv_head_0/bn/betasave_1/RestoreV2:60*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/bn/beta*
_output_shapes
:@
ä
save_1/Assign_61Assign(layerfilter1_newfea_conv_head_0/bn/gammasave_1/RestoreV2:61*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_newfea_conv_head_0/bn/gamma*
_output_shapes
:@
ñ
save_1/Assign_62Assignnlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:62*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
õ
save_1/Assign_63Assignplayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:63*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
â
save_1/Assign_64Assign'layerfilter1_newfea_conv_head_1/bn/betasave_1/RestoreV2:64*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/bn/beta*
_output_shapes
:@
ä
save_1/Assign_65Assign(layerfilter1_newfea_conv_head_1/bn/gammasave_1/RestoreV2:65*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_newfea_conv_head_1/bn/gamma*
_output_shapes
:@
ñ
save_1/Assign_66Assignnlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:66*
T0*
validate_shape(*
use_locking(*
_classw
usloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:@
õ
save_1/Assign_67Assignplayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:67*
T0*
validate_shape(*
use_locking(*
_classy
wuloc:@layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:@
æ
save_1/Assign_68Assign)layerfilter1_self_att_conv_head_0/bn/betasave_1/RestoreV2:68*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/bn/beta*
_output_shapes
:
è
save_1/Assign_69Assign*layerfilter1_self_att_conv_head_0/bn/gammasave_1/RestoreV2:69*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_self_att_conv_head_0/bn/gamma*
_output_shapes
:
ù
save_1/Assign_70Assignrlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:70*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ý
save_1/Assign_71Assigntlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:71*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
æ
save_1/Assign_72Assign)layerfilter1_self_att_conv_head_1/bn/betasave_1/RestoreV2:72*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/bn/beta*
_output_shapes
:
è
save_1/Assign_73Assign*layerfilter1_self_att_conv_head_1/bn/gammasave_1/RestoreV2:73*
T0*
validate_shape(*
use_locking(*=
_class3
1/loc:@layerfilter1_self_att_conv_head_1/bn/gamma*
_output_shapes
:
ù
save_1/Assign_74Assignrlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAveragesave_1/RestoreV2:74*
T0*
validate_shape(*
use_locking(*
_class{
ywloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage*
_output_shapes
:
ý
save_1/Assign_75Assigntlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAveragesave_1/RestoreV2:75*
T0*
validate_shape(*
use_locking(*
_class}
{yloc:@layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage*
_output_shapes
:
´
save_1/restore_shardNoOp^save_1/Assign^save_1/Assign_1^save_1/Assign_2^save_1/Assign_3^save_1/Assign_4^save_1/Assign_5^save_1/Assign_6^save_1/Assign_7^save_1/Assign_8^save_1/Assign_9^save_1/Assign_10^save_1/Assign_11^save_1/Assign_12^save_1/Assign_13^save_1/Assign_14^save_1/Assign_15^save_1/Assign_16^save_1/Assign_17^save_1/Assign_18^save_1/Assign_19^save_1/Assign_20^save_1/Assign_21^save_1/Assign_22^save_1/Assign_23^save_1/Assign_24^save_1/Assign_25^save_1/Assign_26^save_1/Assign_27^save_1/Assign_28^save_1/Assign_29^save_1/Assign_30^save_1/Assign_31^save_1/Assign_32^save_1/Assign_33^save_1/Assign_34^save_1/Assign_35^save_1/Assign_36^save_1/Assign_37^save_1/Assign_38^save_1/Assign_39^save_1/Assign_40^save_1/Assign_41^save_1/Assign_42^save_1/Assign_43^save_1/Assign_44^save_1/Assign_45^save_1/Assign_46^save_1/Assign_47^save_1/Assign_48^save_1/Assign_49^save_1/Assign_50^save_1/Assign_51^save_1/Assign_52^save_1/Assign_53^save_1/Assign_54^save_1/Assign_55^save_1/Assign_56^save_1/Assign_57^save_1/Assign_58^save_1/Assign_59^save_1/Assign_60^save_1/Assign_61^save_1/Assign_62^save_1/Assign_63^save_1/Assign_64^save_1/Assign_65^save_1/Assign_66^save_1/Assign_67^save_1/Assign_68^save_1/Assign_69^save_1/Assign_70^save_1/Assign_71^save_1/Assign_72^save_1/Assign_73^save_1/Assign_74^save_1/Assign_75

save_1/RestoreV2_1/tensor_namesConst"/device:CPU:0*´

valueª
B§
/B
agg/biasesBagg/weightsBgapnet00/biasesBgapnet00/weightsBgapnet01/biasesBgapnet01/weightsBgapnet10/biasesBgapnet10/weightsBgapnet11/biasesBgapnet11/weightsBglobal_expand/biasesBglobal_expand/weightsBlayerfilter0_edgefea_0/biasesBlayerfilter0_edgefea_0/weightsB(layerfilter0_neib_att_conv_head_0/biasesB)layerfilter0_neib_att_conv_head_0/weightsB'layerfilter0_newfea_conv_head_0/weightsB(layerfilter0_self_att_conv_head_0/biasesB)layerfilter0_self_att_conv_head_0/weightsBlayerfilter1_edgefea_0/biasesBlayerfilter1_edgefea_0/weightsBlayerfilter1_edgefea_1/biasesBlayerfilter1_edgefea_1/weightsB(layerfilter1_neib_att_conv_head_0/biasesB)layerfilter1_neib_att_conv_head_0/weightsB(layerfilter1_neib_att_conv_head_1/biasesB)layerfilter1_neib_att_conv_head_1/weightsB'layerfilter1_newfea_conv_head_0/weightsB'layerfilter1_newfea_conv_head_1/weightsB(layerfilter1_self_att_conv_head_0/biasesB)layerfilter1_self_att_conv_head_0/weightsB(layerfilter1_self_att_conv_head_1/biasesB)layerfilter1_self_att_conv_head_1/weightsBseg/conv2/biasesBseg/conv2/bn/betaBseg/conv2/bn/gammaBseg/conv2/bn/pop_meanBseg/conv2/bn/pop_varBseg/conv2/weightsBseg/conv3/biasesBseg/conv3/bn/betaBseg/conv3/bn/gammaBseg/conv3/bn/pop_meanBseg/conv3/bn/pop_varBseg/conv3/weightsBseg/conv5/biasesBseg/conv5/weights*
dtype0*
_output_shapes
:/
×
#save_1/RestoreV2_1/shape_and_slicesConst"/device:CPU:0*q
valuehBf/B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B *
dtype0*
_output_shapes
:/

save_1/RestoreV2_1	RestoreV2save_1/Constsave_1/RestoreV2_1/tensor_names#save_1/RestoreV2_1/shape_and_slices"/device:CPU:0*=
dtypes3
12/*Ò
_output_shapes¿
¼:::::::::::::::::::::::::::::::::::::::::::::::
·
save_1/Assign_76Assign
agg/biasessave_1/RestoreV2_1"/device:CPU:0*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/biases*
_output_shapes	
:
È
save_1/Assign_77Assignagg/weightssave_1/RestoreV2_1:1"/device:CPU:0*
T0*
validate_shape(*
use_locking(*
_class
loc:@agg/weights*(
_output_shapes
:ð
Â
save_1/Assign_78Assigngapnet00/biasessave_1/RestoreV2_1:2"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet00/biases*
_output_shapes
:@
Ð
save_1/Assign_79Assigngapnet00/weightssave_1/RestoreV2_1:3"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet00/weights*&
_output_shapes
:(@
Ã
save_1/Assign_80Assigngapnet01/biasessave_1/RestoreV2_1:4"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet01/biases*
_output_shapes	
:
Ñ
save_1/Assign_81Assigngapnet01/weightssave_1/RestoreV2_1:5"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet01/weights*'
_output_shapes
:@
Ã
save_1/Assign_82Assigngapnet10/biasessave_1/RestoreV2_1:6"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet10/biases*
_output_shapes	
:
Ò
save_1/Assign_83Assigngapnet10/weightssave_1/RestoreV2_1:7"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet10/weights*(
_output_shapes
:
Ã
save_1/Assign_84Assigngapnet11/biasessave_1/RestoreV2_1:8"/device:CPU:0*
T0*
validate_shape(*
use_locking(*"
_class
loc:@gapnet11/biases*
_output_shapes	
:
Ò
save_1/Assign_85Assigngapnet11/weightssave_1/RestoreV2_1:9"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@gapnet11/weights*(
_output_shapes
:
Í
save_1/Assign_86Assignglobal_expand/biasessave_1/RestoreV2_1:10"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@global_expand/biases*
_output_shapes
:
Û
save_1/Assign_87Assignglobal_expand/weightssave_1/RestoreV2_1:11"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@global_expand/weights*&
_output_shapes
:
ß
save_1/Assign_88Assignlayerfilter0_edgefea_0/biasessave_1/RestoreV2_1:12"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter0_edgefea_0/biases*
_output_shapes
: 
í
save_1/Assign_89Assignlayerfilter0_edgefea_0/weightssave_1/RestoreV2_1:13"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter0_edgefea_0/weights*&
_output_shapes
: 
õ
save_1/Assign_90Assign(layerfilter0_neib_att_conv_head_0/biasessave_1/RestoreV2_1:14"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_neib_att_conv_head_0/biases*
_output_shapes
:

save_1/Assign_91Assign)layerfilter0_neib_att_conv_head_0/weightssave_1/RestoreV2_1:15"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_neib_att_conv_head_0/weights*&
_output_shapes
: 
ÿ
save_1/Assign_92Assign'layerfilter0_newfea_conv_head_0/weightssave_1/RestoreV2_1:16"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter0_newfea_conv_head_0/weights*&
_output_shapes
: 
õ
save_1/Assign_93Assign(layerfilter0_self_att_conv_head_0/biasessave_1/RestoreV2_1:17"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter0_self_att_conv_head_0/biases*
_output_shapes
:

save_1/Assign_94Assign)layerfilter0_self_att_conv_head_0/weightssave_1/RestoreV2_1:18"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter0_self_att_conv_head_0/weights*&
_output_shapes
: 
ß
save_1/Assign_95Assignlayerfilter1_edgefea_0/biasessave_1/RestoreV2_1:19"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter1_edgefea_0/biases*
_output_shapes
:@
î
save_1/Assign_96Assignlayerfilter1_edgefea_0/weightssave_1/RestoreV2_1:20"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_0/weights*'
_output_shapes
:@
ß
save_1/Assign_97Assignlayerfilter1_edgefea_1/biasessave_1/RestoreV2_1:21"/device:CPU:0*
T0*
validate_shape(*
use_locking(*0
_class&
$"loc:@layerfilter1_edgefea_1/biases*
_output_shapes
:@
î
save_1/Assign_98Assignlayerfilter1_edgefea_1/weightssave_1/RestoreV2_1:22"/device:CPU:0*
T0*
validate_shape(*
use_locking(*1
_class'
%#loc:@layerfilter1_edgefea_1/weights*'
_output_shapes
:@
õ
save_1/Assign_99Assign(layerfilter1_neib_att_conv_head_0/biasessave_1/RestoreV2_1:23"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_0/biases*
_output_shapes
:

save_1/Assign_100Assign)layerfilter1_neib_att_conv_head_0/weightssave_1/RestoreV2_1:24"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_0/weights*&
_output_shapes
:@
ö
save_1/Assign_101Assign(layerfilter1_neib_att_conv_head_1/biasessave_1/RestoreV2_1:25"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_neib_att_conv_head_1/biases*
_output_shapes
:

save_1/Assign_102Assign)layerfilter1_neib_att_conv_head_1/weightssave_1/RestoreV2_1:26"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_neib_att_conv_head_1/weights*&
_output_shapes
:@

save_1/Assign_103Assign'layerfilter1_newfea_conv_head_0/weightssave_1/RestoreV2_1:27"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_0/weights*'
_output_shapes
:@

save_1/Assign_104Assign'layerfilter1_newfea_conv_head_1/weightssave_1/RestoreV2_1:28"/device:CPU:0*
T0*
validate_shape(*
use_locking(*:
_class0
.,loc:@layerfilter1_newfea_conv_head_1/weights*'
_output_shapes
:@
ö
save_1/Assign_105Assign(layerfilter1_self_att_conv_head_0/biasessave_1/RestoreV2_1:29"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_self_att_conv_head_0/biases*
_output_shapes
:

save_1/Assign_106Assign)layerfilter1_self_att_conv_head_0/weightssave_1/RestoreV2_1:30"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_0/weights*&
_output_shapes
:@
ö
save_1/Assign_107Assign(layerfilter1_self_att_conv_head_1/biasessave_1/RestoreV2_1:31"/device:CPU:0*
T0*
validate_shape(*
use_locking(*;
_class1
/-loc:@layerfilter1_self_att_conv_head_1/biases*
_output_shapes
:

save_1/Assign_108Assign)layerfilter1_self_att_conv_head_1/weightssave_1/RestoreV2_1:32"/device:CPU:0*
T0*
validate_shape(*
use_locking(*<
_class2
0.loc:@layerfilter1_self_att_conv_head_1/weights*&
_output_shapes
:@
Ç
save_1/Assign_109Assignseg/conv2/biasessave_1/RestoreV2_1:33"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv2/biases*
_output_shapes	
:
É
save_1/Assign_110Assignseg/conv2/bn/betasave_1/RestoreV2_1:34"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv2/bn/beta*
_output_shapes	
:
Ë
save_1/Assign_111Assignseg/conv2/bn/gammasave_1/RestoreV2_1:35"/device:CPU:0*
T0*
validate_shape(*
use_locking(*%
_class
loc:@seg/conv2/bn/gamma*
_output_shapes	
:
Ñ
save_1/Assign_112Assignseg/conv2/bn/pop_meansave_1/RestoreV2_1:36"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv2/bn/pop_mean*
_output_shapes	
:
Ï
save_1/Assign_113Assignseg/conv2/bn/pop_varsave_1/RestoreV2_1:37"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv2/bn/pop_var*
_output_shapes	
:
Ö
save_1/Assign_114Assignseg/conv2/weightssave_1/RestoreV2_1:38"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv2/weights*(
_output_shapes
:
Ç
save_1/Assign_115Assignseg/conv3/biasessave_1/RestoreV2_1:39"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv3/biases*
_output_shapes	
:
É
save_1/Assign_116Assignseg/conv3/bn/betasave_1/RestoreV2_1:40"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv3/bn/beta*
_output_shapes	
:
Ë
save_1/Assign_117Assignseg/conv3/bn/gammasave_1/RestoreV2_1:41"/device:CPU:0*
T0*
validate_shape(*
use_locking(*%
_class
loc:@seg/conv3/bn/gamma*
_output_shapes	
:
Ñ
save_1/Assign_118Assignseg/conv3/bn/pop_meansave_1/RestoreV2_1:42"/device:CPU:0*
T0*
validate_shape(*
use_locking(*(
_class
loc:@seg/conv3/bn/pop_mean*
_output_shapes	
:
Ï
save_1/Assign_119Assignseg/conv3/bn/pop_varsave_1/RestoreV2_1:43"/device:CPU:0*
T0*
validate_shape(*
use_locking(*'
_class
loc:@seg/conv3/bn/pop_var*
_output_shapes	
:
Ö
save_1/Assign_120Assignseg/conv3/weightssave_1/RestoreV2_1:44"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv3/weights*(
_output_shapes
:
Æ
save_1/Assign_121Assignseg/conv5/biasessave_1/RestoreV2_1:45"/device:CPU:0*
T0*
validate_shape(*
use_locking(*#
_class
loc:@seg/conv5/biases*
_output_shapes
:
Õ
save_1/Assign_122Assignseg/conv5/weightssave_1/RestoreV2_1:46"/device:CPU:0*
T0*
validate_shape(*
use_locking(*$
_class
loc:@seg/conv5/weights*'
_output_shapes
:
Á
save_1/restore_shard_1NoOp^save_1/Assign_76^save_1/Assign_77^save_1/Assign_78^save_1/Assign_79^save_1/Assign_80^save_1/Assign_81^save_1/Assign_82^save_1/Assign_83^save_1/Assign_84^save_1/Assign_85^save_1/Assign_86^save_1/Assign_87^save_1/Assign_88^save_1/Assign_89^save_1/Assign_90^save_1/Assign_91^save_1/Assign_92^save_1/Assign_93^save_1/Assign_94^save_1/Assign_95^save_1/Assign_96^save_1/Assign_97^save_1/Assign_98^save_1/Assign_99^save_1/Assign_100^save_1/Assign_101^save_1/Assign_102^save_1/Assign_103^save_1/Assign_104^save_1/Assign_105^save_1/Assign_106^save_1/Assign_107^save_1/Assign_108^save_1/Assign_109^save_1/Assign_110^save_1/Assign_111^save_1/Assign_112^save_1/Assign_113^save_1/Assign_114^save_1/Assign_115^save_1/Assign_116^save_1/Assign_117^save_1/Assign_118^save_1/Assign_119^save_1/Assign_120^save_1/Assign_121^save_1/Assign_122"/device:CPU:0
6
save_1/restore_all/NoOpNoOp^save_1/restore_shard
I
save_1/restore_all/NoOp_1NoOp^save_1/restore_shard_1"/device:CPU:0
P
save_1/restore_allNoOp^save_1/restore_all/NoOp^save_1/restore_all/NoOp_1"B
save_1/Const:0save_1/Identity:0save_1/restore_all (5 @F8"¢Þ
	variablesÞÞ
H

Variable:0Variable/AssignVariable/read:02Variable/initial_value:0
Ñ
)layerfilter0_newfea_conv_head_0/weights:0.layerfilter0_newfea_conv_head_0/weights/Assign.layerfilter0_newfea_conv_head_0/weights/read:02Dlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform:0
·
)layerfilter0_newfea_conv_head_0/bn/beta:0.layerfilter0_newfea_conv_head_0/bn/beta/Assign.layerfilter0_newfea_conv_head_0/bn/beta/read:02*layerfilter0_newfea_conv_head_0/bn/Const:0
¼
*layerfilter0_newfea_conv_head_0/bn/gamma:0/layerfilter0_newfea_conv_head_0/bn/gamma/Assign/layerfilter0_newfea_conv_head_0/bn/gamma/read:02,layerfilter0_newfea_conv_head_0/bn/Const_1:0
å
playerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
í
rlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0wlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignwlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
­
 layerfilter0_edgefea_0/weights:0%layerfilter0_edgefea_0/weights/Assign%layerfilter0_edgefea_0/weights/read:02;layerfilter0_edgefea_0/weights/Initializer/random_uniform:0
 
layerfilter0_edgefea_0/biases:0$layerfilter0_edgefea_0/biases/Assign$layerfilter0_edgefea_0/biases/read:021layerfilter0_edgefea_0/biases/Initializer/Const:0

 layerfilter0_edgefea_0/bn/beta:0%layerfilter0_edgefea_0/bn/beta/Assign%layerfilter0_edgefea_0/bn/beta/read:02!layerfilter0_edgefea_0/bn/Const:0

!layerfilter0_edgefea_0/bn/gamma:0&layerfilter0_edgefea_0/bn/gamma/Assign&layerfilter0_edgefea_0/bn/gamma/read:02#layerfilter0_edgefea_0/bn/Const_1:0

^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0clayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignclayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02playerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
¤
`layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0elayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignelayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02rlayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Ù
+layerfilter0_self_att_conv_head_0/weights:00layerfilter0_self_att_conv_head_0/weights/Assign0layerfilter0_self_att_conv_head_0/weights/read:02Flayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter0_self_att_conv_head_0/biases:0/layerfilter0_self_att_conv_head_0/biases/Assign/layerfilter0_self_att_conv_head_0/biases/read:02<layerfilter0_self_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter0_self_att_conv_head_0/bn/beta:00layerfilter0_self_att_conv_head_0/bn/beta/Assign0layerfilter0_self_att_conv_head_0/bn/beta/read:02,layerfilter0_self_att_conv_head_0/bn/Const:0
Ä
,layerfilter0_self_att_conv_head_0/bn/gamma:01layerfilter0_self_att_conv_head_0/bn/gamma/Assign1layerfilter0_self_att_conv_head_0/bn/gamma/read:02.layerfilter0_self_att_conv_head_0/bn/Const_1:0
õ
tlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Ù
+layerfilter0_neib_att_conv_head_0/weights:00layerfilter0_neib_att_conv_head_0/weights/Assign0layerfilter0_neib_att_conv_head_0/weights/read:02Flayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter0_neib_att_conv_head_0/biases:0/layerfilter0_neib_att_conv_head_0/biases/Assign/layerfilter0_neib_att_conv_head_0/biases/read:02<layerfilter0_neib_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter0_neib_att_conv_head_0/bn/beta:00layerfilter0_neib_att_conv_head_0/bn/beta/Assign0layerfilter0_neib_att_conv_head_0/bn/beta/read:02,layerfilter0_neib_att_conv_head_0/bn/Const:0
Ä
,layerfilter0_neib_att_conv_head_0/bn/gamma:01layerfilter0_neib_att_conv_head_0/bn/gamma/Assign1layerfilter0_neib_att_conv_head_0/bn/gamma/read:02.layerfilter0_neib_att_conv_head_0/bn/Const_1:0
õ
tlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
d
BiasAdd/biases:0BiasAdd/biases/AssignBiasAdd/biases/read:02"BiasAdd/biases/Initializer/zeros:0
u
gapnet00/weights:0gapnet00/weights/Assigngapnet00/weights/read:02-gapnet00/weights/Initializer/random_uniform:0
h
gapnet00/biases:0gapnet00/biases/Assigngapnet00/biases/read:02#gapnet00/biases/Initializer/Const:0
[
gapnet00/bn/beta:0gapnet00/bn/beta/Assigngapnet00/bn/beta/read:02gapnet00/bn/Const:0
`
gapnet00/bn/gamma:0gapnet00/bn/gamma/Assigngapnet00/bn/gamma/read:02gapnet00/bn/Const_1:0
¬
Bgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage:0Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/AssignGgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:02Tgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
Dgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage:0Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignIgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Vgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
u
gapnet01/weights:0gapnet01/weights/Assigngapnet01/weights/read:02-gapnet01/weights/Initializer/random_uniform:0
h
gapnet01/biases:0gapnet01/biases/Assigngapnet01/biases/read:02#gapnet01/biases/Initializer/Const:0
[
gapnet01/bn/beta:0gapnet01/bn/beta/Assigngapnet01/bn/beta/read:02gapnet01/bn/Const:0
`
gapnet01/bn/gamma:0gapnet01/bn/gamma/Assigngapnet01/bn/gamma/read:02gapnet01/bn/Const_1:0
¬
Bgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage:0Ggapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/AssignGgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/read:02Tgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
Dgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage:0Igapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignIgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Vgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Ñ
)layerfilter1_newfea_conv_head_0/weights:0.layerfilter1_newfea_conv_head_0/weights/Assign.layerfilter1_newfea_conv_head_0/weights/read:02Dlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform:0
·
)layerfilter1_newfea_conv_head_0/bn/beta:0.layerfilter1_newfea_conv_head_0/bn/beta/Assign.layerfilter1_newfea_conv_head_0/bn/beta/read:02*layerfilter1_newfea_conv_head_0/bn/Const:0
¼
*layerfilter1_newfea_conv_head_0/bn/gamma:0/layerfilter1_newfea_conv_head_0/bn/gamma/Assign/layerfilter1_newfea_conv_head_0/bn/gamma/read:02,layerfilter1_newfea_conv_head_0/bn/Const_1:0
å
playerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
í
rlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0wlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignwlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
­
 layerfilter1_edgefea_0/weights:0%layerfilter1_edgefea_0/weights/Assign%layerfilter1_edgefea_0/weights/read:02;layerfilter1_edgefea_0/weights/Initializer/random_uniform:0
 
layerfilter1_edgefea_0/biases:0$layerfilter1_edgefea_0/biases/Assign$layerfilter1_edgefea_0/biases/read:021layerfilter1_edgefea_0/biases/Initializer/Const:0

 layerfilter1_edgefea_0/bn/beta:0%layerfilter1_edgefea_0/bn/beta/Assign%layerfilter1_edgefea_0/bn/beta/read:02!layerfilter1_edgefea_0/bn/Const:0

!layerfilter1_edgefea_0/bn/gamma:0&layerfilter1_edgefea_0/bn/gamma/Assign&layerfilter1_edgefea_0/bn/gamma/read:02#layerfilter1_edgefea_0/bn/Const_1:0

^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0clayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignclayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02playerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
¤
`layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0elayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignelayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02rlayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Ù
+layerfilter1_self_att_conv_head_0/weights:00layerfilter1_self_att_conv_head_0/weights/Assign0layerfilter1_self_att_conv_head_0/weights/read:02Flayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter1_self_att_conv_head_0/biases:0/layerfilter1_self_att_conv_head_0/biases/Assign/layerfilter1_self_att_conv_head_0/biases/read:02<layerfilter1_self_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter1_self_att_conv_head_0/bn/beta:00layerfilter1_self_att_conv_head_0/bn/beta/Assign0layerfilter1_self_att_conv_head_0/bn/beta/read:02,layerfilter1_self_att_conv_head_0/bn/Const:0
Ä
,layerfilter1_self_att_conv_head_0/bn/gamma:01layerfilter1_self_att_conv_head_0/bn/gamma/Assign1layerfilter1_self_att_conv_head_0/bn/gamma/read:02.layerfilter1_self_att_conv_head_0/bn/Const_1:0
õ
tlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Ù
+layerfilter1_neib_att_conv_head_0/weights:00layerfilter1_neib_att_conv_head_0/weights/Assign0layerfilter1_neib_att_conv_head_0/weights/read:02Flayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter1_neib_att_conv_head_0/biases:0/layerfilter1_neib_att_conv_head_0/biases/Assign/layerfilter1_neib_att_conv_head_0/biases/read:02<layerfilter1_neib_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter1_neib_att_conv_head_0/bn/beta:00layerfilter1_neib_att_conv_head_0/bn/beta/Assign0layerfilter1_neib_att_conv_head_0/bn/beta/read:02,layerfilter1_neib_att_conv_head_0/bn/Const:0
Ä
,layerfilter1_neib_att_conv_head_0/bn/gamma:01layerfilter1_neib_att_conv_head_0/bn/gamma/Assign1layerfilter1_neib_att_conv_head_0/bn/gamma/read:02.layerfilter1_neib_att_conv_head_0/bn/Const_1:0
õ
tlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
l
BiasAdd_1/biases:0BiasAdd_1/biases/AssignBiasAdd_1/biases/read:02$BiasAdd_1/biases/Initializer/zeros:0
Ñ
)layerfilter1_newfea_conv_head_1/weights:0.layerfilter1_newfea_conv_head_1/weights/Assign.layerfilter1_newfea_conv_head_1/weights/read:02Dlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform:0
·
)layerfilter1_newfea_conv_head_1/bn/beta:0.layerfilter1_newfea_conv_head_1/bn/beta/Assign.layerfilter1_newfea_conv_head_1/bn/beta/read:02*layerfilter1_newfea_conv_head_1/bn/Const:0
¼
*layerfilter1_newfea_conv_head_1/bn/gamma:0/layerfilter1_newfea_conv_head_1/bn/gamma/Assign/layerfilter1_newfea_conv_head_1/bn/gamma/read:02,layerfilter1_newfea_conv_head_1/bn/Const_1:0
å
playerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0ulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Assignulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
í
rlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0wlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignwlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
­
 layerfilter1_edgefea_1/weights:0%layerfilter1_edgefea_1/weights/Assign%layerfilter1_edgefea_1/weights/read:02;layerfilter1_edgefea_1/weights/Initializer/random_uniform:0
 
layerfilter1_edgefea_1/biases:0$layerfilter1_edgefea_1/biases/Assign$layerfilter1_edgefea_1/biases/read:021layerfilter1_edgefea_1/biases/Initializer/Const:0

 layerfilter1_edgefea_1/bn/beta:0%layerfilter1_edgefea_1/bn/beta/Assign%layerfilter1_edgefea_1/bn/beta/read:02!layerfilter1_edgefea_1/bn/Const:0

!layerfilter1_edgefea_1/bn/gamma:0&layerfilter1_edgefea_1/bn/gamma/Assign&layerfilter1_edgefea_1/bn/gamma/read:02#layerfilter1_edgefea_1/bn/Const_1:0

^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage:0clayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Assignclayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/read:02playerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
¤
`layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0elayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Assignelayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02rlayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Ù
+layerfilter1_self_att_conv_head_1/weights:00layerfilter1_self_att_conv_head_1/weights/Assign0layerfilter1_self_att_conv_head_1/weights/read:02Flayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform:0
Ì
*layerfilter1_self_att_conv_head_1/biases:0/layerfilter1_self_att_conv_head_1/biases/Assign/layerfilter1_self_att_conv_head_1/biases/read:02<layerfilter1_self_att_conv_head_1/biases/Initializer/Const:0
¿
+layerfilter1_self_att_conv_head_1/bn/beta:00layerfilter1_self_att_conv_head_1/bn/beta/Assign0layerfilter1_self_att_conv_head_1/bn/beta/read:02,layerfilter1_self_att_conv_head_1/bn/Const:0
Ä
,layerfilter1_self_att_conv_head_1/bn/gamma:01layerfilter1_self_att_conv_head_1/bn/gamma/Assign1layerfilter1_self_att_conv_head_1/bn/gamma/read:02.layerfilter1_self_att_conv_head_1/bn/Const_1:0
õ
tlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
Ù
+layerfilter1_neib_att_conv_head_1/weights:00layerfilter1_neib_att_conv_head_1/weights/Assign0layerfilter1_neib_att_conv_head_1/weights/read:02Flayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform:0
Ì
*layerfilter1_neib_att_conv_head_1/biases:0/layerfilter1_neib_att_conv_head_1/biases/Assign/layerfilter1_neib_att_conv_head_1/biases/read:02<layerfilter1_neib_att_conv_head_1/biases/Initializer/Const:0
¿
+layerfilter1_neib_att_conv_head_1/bn/beta:00layerfilter1_neib_att_conv_head_1/bn/beta/Assign0layerfilter1_neib_att_conv_head_1/bn/beta/read:02,layerfilter1_neib_att_conv_head_1/bn/Const:0
Ä
,layerfilter1_neib_att_conv_head_1/bn/gamma:01layerfilter1_neib_att_conv_head_1/bn/gamma/Assign1layerfilter1_neib_att_conv_head_1/bn/gamma/read:02.layerfilter1_neib_att_conv_head_1/bn/Const_1:0
õ
tlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0ylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Assignylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:02layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
ý
vlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0{layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign{layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
l
BiasAdd_2/biases:0BiasAdd_2/biases/AssignBiasAdd_2/biases/read:02$BiasAdd_2/biases/Initializer/zeros:0
u
gapnet10/weights:0gapnet10/weights/Assigngapnet10/weights/read:02-gapnet10/weights/Initializer/random_uniform:0
h
gapnet10/biases:0gapnet10/biases/Assigngapnet10/biases/read:02#gapnet10/biases/Initializer/Const:0
[
gapnet10/bn/beta:0gapnet10/bn/beta/Assigngapnet10/bn/beta/read:02gapnet10/bn/Const:0
`
gapnet10/bn/gamma:0gapnet10/bn/gamma/Assigngapnet10/bn/gamma/read:02gapnet10/bn/Const_1:0
¬
Bgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage:0Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/AssignGgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:02Tgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
Dgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage:0Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignIgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Vgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
u
gapnet11/weights:0gapnet11/weights/Assigngapnet11/weights/read:02-gapnet11/weights/Initializer/random_uniform:0
h
gapnet11/biases:0gapnet11/biases/Assigngapnet11/biases/read:02#gapnet11/biases/Initializer/Const:0
[
gapnet11/bn/beta:0gapnet11/bn/beta/Assigngapnet11/bn/beta/read:02gapnet11/bn/Const:0
`
gapnet11/bn/gamma:0gapnet11/bn/gamma/Assigngapnet11/bn/gamma/read:02gapnet11/bn/Const_1:0
¬
Bgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage:0Ggapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/AssignGgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/read:02Tgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
´
Dgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage:0Igapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignIgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Vgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0

global_expand/weights:0global_expand/weights/Assignglobal_expand/weights/read:022global_expand/weights/Initializer/random_uniform:0
|
global_expand/biases:0global_expand/biases/Assignglobal_expand/biases/read:02(global_expand/biases/Initializer/Const:0
o
global_expand/bn/beta:0global_expand/bn/beta/Assignglobal_expand/bn/beta/read:02global_expand/bn/Const:0
t
global_expand/bn/gamma:0global_expand/bn/gamma/Assignglobal_expand/bn/gamma/read:02global_expand/bn/Const_1:0
Ô
Lglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage:0Qglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/AssignQglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/read:02^global_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0
Ü
Nglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage:0Sglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/AssignSglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02`global_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
a
agg/weights:0agg/weights/Assignagg/weights/read:02(agg/weights/Initializer/random_uniform:0
T
agg/biases:0agg/biases/Assignagg/biases/read:02agg/biases/Initializer/Const:0
G
agg/bn/beta:0agg/bn/beta/Assignagg/bn/beta/read:02agg/bn/Const:0
L
agg/bn/gamma:0agg/bn/gamma/Assignagg/bn/gamma/read:02agg/bn/Const_1:0

8agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage:0=agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Assign=agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/read:02Jagg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/Initializer/zeros:0

:agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage:0?agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Assign?agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/read:02Lagg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/Initializer/zeros:0
y
seg/conv2/weights:0seg/conv2/weights/Assignseg/conv2/weights/read:02.seg/conv2/weights/Initializer/random_uniform:0
l
seg/conv2/biases:0seg/conv2/biases/Assignseg/conv2/biases/read:02$seg/conv2/biases/Initializer/Const:0
p
seg/conv2/bn/beta:0seg/conv2/bn/beta/Assignseg/conv2/bn/beta/read:02%seg/conv2/bn/beta/Initializer/zeros:0
s
seg/conv2/bn/gamma:0seg/conv2/bn/gamma/Assignseg/conv2/bn/gamma/read:02%seg/conv2/bn/gamma/Initializer/ones:0

seg/conv2/bn/pop_mean:0seg/conv2/bn/pop_mean/Assignseg/conv2/bn/pop_mean/read:02)seg/conv2/bn/pop_mean/Initializer/zeros:0
{
seg/conv2/bn/pop_var:0seg/conv2/bn/pop_var/Assignseg/conv2/bn/pop_var/read:02'seg/conv2/bn/pop_var/Initializer/ones:0
y
seg/conv3/weights:0seg/conv3/weights/Assignseg/conv3/weights/read:02.seg/conv3/weights/Initializer/random_uniform:0
l
seg/conv3/biases:0seg/conv3/biases/Assignseg/conv3/biases/read:02$seg/conv3/biases/Initializer/Const:0
p
seg/conv3/bn/beta:0seg/conv3/bn/beta/Assignseg/conv3/bn/beta/read:02%seg/conv3/bn/beta/Initializer/zeros:0
s
seg/conv3/bn/gamma:0seg/conv3/bn/gamma/Assignseg/conv3/bn/gamma/read:02%seg/conv3/bn/gamma/Initializer/ones:0

seg/conv3/bn/pop_mean:0seg/conv3/bn/pop_mean/Assignseg/conv3/bn/pop_mean/read:02)seg/conv3/bn/pop_mean/Initializer/zeros:0
{
seg/conv3/bn/pop_var:0seg/conv3/bn/pop_var/Assignseg/conv3/bn/pop_var/read:02'seg/conv3/bn/pop_var/Initializer/ones:0
y
seg/conv5/weights:0seg/conv5/weights/Assignseg/conv5/weights/read:02.seg/conv5/weights/Initializer/random_uniform:0
l
seg/conv5/biases:0seg/conv5/biases/Assignseg/conv5/biases/read:02$seg/conv5/biases/Initializer/Const:0"Òb
trainable_variablesºb·b
Ñ
)layerfilter0_newfea_conv_head_0/weights:0.layerfilter0_newfea_conv_head_0/weights/Assign.layerfilter0_newfea_conv_head_0/weights/read:02Dlayerfilter0_newfea_conv_head_0/weights/Initializer/random_uniform:0
·
)layerfilter0_newfea_conv_head_0/bn/beta:0.layerfilter0_newfea_conv_head_0/bn/beta/Assign.layerfilter0_newfea_conv_head_0/bn/beta/read:02*layerfilter0_newfea_conv_head_0/bn/Const:0
¼
*layerfilter0_newfea_conv_head_0/bn/gamma:0/layerfilter0_newfea_conv_head_0/bn/gamma/Assign/layerfilter0_newfea_conv_head_0/bn/gamma/read:02,layerfilter0_newfea_conv_head_0/bn/Const_1:0
­
 layerfilter0_edgefea_0/weights:0%layerfilter0_edgefea_0/weights/Assign%layerfilter0_edgefea_0/weights/read:02;layerfilter0_edgefea_0/weights/Initializer/random_uniform:0
 
layerfilter0_edgefea_0/biases:0$layerfilter0_edgefea_0/biases/Assign$layerfilter0_edgefea_0/biases/read:021layerfilter0_edgefea_0/biases/Initializer/Const:0

 layerfilter0_edgefea_0/bn/beta:0%layerfilter0_edgefea_0/bn/beta/Assign%layerfilter0_edgefea_0/bn/beta/read:02!layerfilter0_edgefea_0/bn/Const:0

!layerfilter0_edgefea_0/bn/gamma:0&layerfilter0_edgefea_0/bn/gamma/Assign&layerfilter0_edgefea_0/bn/gamma/read:02#layerfilter0_edgefea_0/bn/Const_1:0
Ù
+layerfilter0_self_att_conv_head_0/weights:00layerfilter0_self_att_conv_head_0/weights/Assign0layerfilter0_self_att_conv_head_0/weights/read:02Flayerfilter0_self_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter0_self_att_conv_head_0/biases:0/layerfilter0_self_att_conv_head_0/biases/Assign/layerfilter0_self_att_conv_head_0/biases/read:02<layerfilter0_self_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter0_self_att_conv_head_0/bn/beta:00layerfilter0_self_att_conv_head_0/bn/beta/Assign0layerfilter0_self_att_conv_head_0/bn/beta/read:02,layerfilter0_self_att_conv_head_0/bn/Const:0
Ä
,layerfilter0_self_att_conv_head_0/bn/gamma:01layerfilter0_self_att_conv_head_0/bn/gamma/Assign1layerfilter0_self_att_conv_head_0/bn/gamma/read:02.layerfilter0_self_att_conv_head_0/bn/Const_1:0
Ù
+layerfilter0_neib_att_conv_head_0/weights:00layerfilter0_neib_att_conv_head_0/weights/Assign0layerfilter0_neib_att_conv_head_0/weights/read:02Flayerfilter0_neib_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter0_neib_att_conv_head_0/biases:0/layerfilter0_neib_att_conv_head_0/biases/Assign/layerfilter0_neib_att_conv_head_0/biases/read:02<layerfilter0_neib_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter0_neib_att_conv_head_0/bn/beta:00layerfilter0_neib_att_conv_head_0/bn/beta/Assign0layerfilter0_neib_att_conv_head_0/bn/beta/read:02,layerfilter0_neib_att_conv_head_0/bn/Const:0
Ä
,layerfilter0_neib_att_conv_head_0/bn/gamma:01layerfilter0_neib_att_conv_head_0/bn/gamma/Assign1layerfilter0_neib_att_conv_head_0/bn/gamma/read:02.layerfilter0_neib_att_conv_head_0/bn/Const_1:0
d
BiasAdd/biases:0BiasAdd/biases/AssignBiasAdd/biases/read:02"BiasAdd/biases/Initializer/zeros:0
u
gapnet00/weights:0gapnet00/weights/Assigngapnet00/weights/read:02-gapnet00/weights/Initializer/random_uniform:0
h
gapnet00/biases:0gapnet00/biases/Assigngapnet00/biases/read:02#gapnet00/biases/Initializer/Const:0
[
gapnet00/bn/beta:0gapnet00/bn/beta/Assigngapnet00/bn/beta/read:02gapnet00/bn/Const:0
`
gapnet00/bn/gamma:0gapnet00/bn/gamma/Assigngapnet00/bn/gamma/read:02gapnet00/bn/Const_1:0
u
gapnet01/weights:0gapnet01/weights/Assigngapnet01/weights/read:02-gapnet01/weights/Initializer/random_uniform:0
h
gapnet01/biases:0gapnet01/biases/Assigngapnet01/biases/read:02#gapnet01/biases/Initializer/Const:0
[
gapnet01/bn/beta:0gapnet01/bn/beta/Assigngapnet01/bn/beta/read:02gapnet01/bn/Const:0
`
gapnet01/bn/gamma:0gapnet01/bn/gamma/Assigngapnet01/bn/gamma/read:02gapnet01/bn/Const_1:0
Ñ
)layerfilter1_newfea_conv_head_0/weights:0.layerfilter1_newfea_conv_head_0/weights/Assign.layerfilter1_newfea_conv_head_0/weights/read:02Dlayerfilter1_newfea_conv_head_0/weights/Initializer/random_uniform:0
·
)layerfilter1_newfea_conv_head_0/bn/beta:0.layerfilter1_newfea_conv_head_0/bn/beta/Assign.layerfilter1_newfea_conv_head_0/bn/beta/read:02*layerfilter1_newfea_conv_head_0/bn/Const:0
¼
*layerfilter1_newfea_conv_head_0/bn/gamma:0/layerfilter1_newfea_conv_head_0/bn/gamma/Assign/layerfilter1_newfea_conv_head_0/bn/gamma/read:02,layerfilter1_newfea_conv_head_0/bn/Const_1:0
­
 layerfilter1_edgefea_0/weights:0%layerfilter1_edgefea_0/weights/Assign%layerfilter1_edgefea_0/weights/read:02;layerfilter1_edgefea_0/weights/Initializer/random_uniform:0
 
layerfilter1_edgefea_0/biases:0$layerfilter1_edgefea_0/biases/Assign$layerfilter1_edgefea_0/biases/read:021layerfilter1_edgefea_0/biases/Initializer/Const:0

 layerfilter1_edgefea_0/bn/beta:0%layerfilter1_edgefea_0/bn/beta/Assign%layerfilter1_edgefea_0/bn/beta/read:02!layerfilter1_edgefea_0/bn/Const:0

!layerfilter1_edgefea_0/bn/gamma:0&layerfilter1_edgefea_0/bn/gamma/Assign&layerfilter1_edgefea_0/bn/gamma/read:02#layerfilter1_edgefea_0/bn/Const_1:0
Ù
+layerfilter1_self_att_conv_head_0/weights:00layerfilter1_self_att_conv_head_0/weights/Assign0layerfilter1_self_att_conv_head_0/weights/read:02Flayerfilter1_self_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter1_self_att_conv_head_0/biases:0/layerfilter1_self_att_conv_head_0/biases/Assign/layerfilter1_self_att_conv_head_0/biases/read:02<layerfilter1_self_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter1_self_att_conv_head_0/bn/beta:00layerfilter1_self_att_conv_head_0/bn/beta/Assign0layerfilter1_self_att_conv_head_0/bn/beta/read:02,layerfilter1_self_att_conv_head_0/bn/Const:0
Ä
,layerfilter1_self_att_conv_head_0/bn/gamma:01layerfilter1_self_att_conv_head_0/bn/gamma/Assign1layerfilter1_self_att_conv_head_0/bn/gamma/read:02.layerfilter1_self_att_conv_head_0/bn/Const_1:0
Ù
+layerfilter1_neib_att_conv_head_0/weights:00layerfilter1_neib_att_conv_head_0/weights/Assign0layerfilter1_neib_att_conv_head_0/weights/read:02Flayerfilter1_neib_att_conv_head_0/weights/Initializer/random_uniform:0
Ì
*layerfilter1_neib_att_conv_head_0/biases:0/layerfilter1_neib_att_conv_head_0/biases/Assign/layerfilter1_neib_att_conv_head_0/biases/read:02<layerfilter1_neib_att_conv_head_0/biases/Initializer/Const:0
¿
+layerfilter1_neib_att_conv_head_0/bn/beta:00layerfilter1_neib_att_conv_head_0/bn/beta/Assign0layerfilter1_neib_att_conv_head_0/bn/beta/read:02,layerfilter1_neib_att_conv_head_0/bn/Const:0
Ä
,layerfilter1_neib_att_conv_head_0/bn/gamma:01layerfilter1_neib_att_conv_head_0/bn/gamma/Assign1layerfilter1_neib_att_conv_head_0/bn/gamma/read:02.layerfilter1_neib_att_conv_head_0/bn/Const_1:0
l
BiasAdd_1/biases:0BiasAdd_1/biases/AssignBiasAdd_1/biases/read:02$BiasAdd_1/biases/Initializer/zeros:0
Ñ
)layerfilter1_newfea_conv_head_1/weights:0.layerfilter1_newfea_conv_head_1/weights/Assign.layerfilter1_newfea_conv_head_1/weights/read:02Dlayerfilter1_newfea_conv_head_1/weights/Initializer/random_uniform:0
·
)layerfilter1_newfea_conv_head_1/bn/beta:0.layerfilter1_newfea_conv_head_1/bn/beta/Assign.layerfilter1_newfea_conv_head_1/bn/beta/read:02*layerfilter1_newfea_conv_head_1/bn/Const:0
¼
*layerfilter1_newfea_conv_head_1/bn/gamma:0/layerfilter1_newfea_conv_head_1/bn/gamma/Assign/layerfilter1_newfea_conv_head_1/bn/gamma/read:02,layerfilter1_newfea_conv_head_1/bn/Const_1:0
­
 layerfilter1_edgefea_1/weights:0%layerfilter1_edgefea_1/weights/Assign%layerfilter1_edgefea_1/weights/read:02;layerfilter1_edgefea_1/weights/Initializer/random_uniform:0
 
layerfilter1_edgefea_1/biases:0$layerfilter1_edgefea_1/biases/Assign$layerfilter1_edgefea_1/biases/read:021layerfilter1_edgefea_1/biases/Initializer/Const:0

 layerfilter1_edgefea_1/bn/beta:0%layerfilter1_edgefea_1/bn/beta/Assign%layerfilter1_edgefea_1/bn/beta/read:02!layerfilter1_edgefea_1/bn/Const:0

!layerfilter1_edgefea_1/bn/gamma:0&layerfilter1_edgefea_1/bn/gamma/Assign&layerfilter1_edgefea_1/bn/gamma/read:02#layerfilter1_edgefea_1/bn/Const_1:0
Ù
+layerfilter1_self_att_conv_head_1/weights:00layerfilter1_self_att_conv_head_1/weights/Assign0layerfilter1_self_att_conv_head_1/weights/read:02Flayerfilter1_self_att_conv_head_1/weights/Initializer/random_uniform:0
Ì
*layerfilter1_self_att_conv_head_1/biases:0/layerfilter1_self_att_conv_head_1/biases/Assign/layerfilter1_self_att_conv_head_1/biases/read:02<layerfilter1_self_att_conv_head_1/biases/Initializer/Const:0
¿
+layerfilter1_self_att_conv_head_1/bn/beta:00layerfilter1_self_att_conv_head_1/bn/beta/Assign0layerfilter1_self_att_conv_head_1/bn/beta/read:02,layerfilter1_self_att_conv_head_1/bn/Const:0
Ä
,layerfilter1_self_att_conv_head_1/bn/gamma:01layerfilter1_self_att_conv_head_1/bn/gamma/Assign1layerfilter1_self_att_conv_head_1/bn/gamma/read:02.layerfilter1_self_att_conv_head_1/bn/Const_1:0
Ù
+layerfilter1_neib_att_conv_head_1/weights:00layerfilter1_neib_att_conv_head_1/weights/Assign0layerfilter1_neib_att_conv_head_1/weights/read:02Flayerfilter1_neib_att_conv_head_1/weights/Initializer/random_uniform:0
Ì
*layerfilter1_neib_att_conv_head_1/biases:0/layerfilter1_neib_att_conv_head_1/biases/Assign/layerfilter1_neib_att_conv_head_1/biases/read:02<layerfilter1_neib_att_conv_head_1/biases/Initializer/Const:0
¿
+layerfilter1_neib_att_conv_head_1/bn/beta:00layerfilter1_neib_att_conv_head_1/bn/beta/Assign0layerfilter1_neib_att_conv_head_1/bn/beta/read:02,layerfilter1_neib_att_conv_head_1/bn/Const:0
Ä
,layerfilter1_neib_att_conv_head_1/bn/gamma:01layerfilter1_neib_att_conv_head_1/bn/gamma/Assign1layerfilter1_neib_att_conv_head_1/bn/gamma/read:02.layerfilter1_neib_att_conv_head_1/bn/Const_1:0
l
BiasAdd_2/biases:0BiasAdd_2/biases/AssignBiasAdd_2/biases/read:02$BiasAdd_2/biases/Initializer/zeros:0
u
gapnet10/weights:0gapnet10/weights/Assigngapnet10/weights/read:02-gapnet10/weights/Initializer/random_uniform:0
h
gapnet10/biases:0gapnet10/biases/Assigngapnet10/biases/read:02#gapnet10/biases/Initializer/Const:0
[
gapnet10/bn/beta:0gapnet10/bn/beta/Assigngapnet10/bn/beta/read:02gapnet10/bn/Const:0
`
gapnet10/bn/gamma:0gapnet10/bn/gamma/Assigngapnet10/bn/gamma/read:02gapnet10/bn/Const_1:0
u
gapnet11/weights:0gapnet11/weights/Assigngapnet11/weights/read:02-gapnet11/weights/Initializer/random_uniform:0
h
gapnet11/biases:0gapnet11/biases/Assigngapnet11/biases/read:02#gapnet11/biases/Initializer/Const:0
[
gapnet11/bn/beta:0gapnet11/bn/beta/Assigngapnet11/bn/beta/read:02gapnet11/bn/Const:0
`
gapnet11/bn/gamma:0gapnet11/bn/gamma/Assigngapnet11/bn/gamma/read:02gapnet11/bn/Const_1:0

global_expand/weights:0global_expand/weights/Assignglobal_expand/weights/read:022global_expand/weights/Initializer/random_uniform:0
|
global_expand/biases:0global_expand/biases/Assignglobal_expand/biases/read:02(global_expand/biases/Initializer/Const:0
o
global_expand/bn/beta:0global_expand/bn/beta/Assignglobal_expand/bn/beta/read:02global_expand/bn/Const:0
t
global_expand/bn/gamma:0global_expand/bn/gamma/Assignglobal_expand/bn/gamma/read:02global_expand/bn/Const_1:0
a
agg/weights:0agg/weights/Assignagg/weights/read:02(agg/weights/Initializer/random_uniform:0
T
agg/biases:0agg/biases/Assignagg/biases/read:02agg/biases/Initializer/Const:0
G
agg/bn/beta:0agg/bn/beta/Assignagg/bn/beta/read:02agg/bn/Const:0
L
agg/bn/gamma:0agg/bn/gamma/Assignagg/bn/gamma/read:02agg/bn/Const_1:0
y
seg/conv2/weights:0seg/conv2/weights/Assignseg/conv2/weights/read:02.seg/conv2/weights/Initializer/random_uniform:0
l
seg/conv2/biases:0seg/conv2/biases/Assignseg/conv2/biases/read:02$seg/conv2/biases/Initializer/Const:0
p
seg/conv2/bn/beta:0seg/conv2/bn/beta/Assignseg/conv2/bn/beta/read:02%seg/conv2/bn/beta/Initializer/zeros:0
s
seg/conv2/bn/gamma:0seg/conv2/bn/gamma/Assignseg/conv2/bn/gamma/read:02%seg/conv2/bn/gamma/Initializer/ones:0
y
seg/conv3/weights:0seg/conv3/weights/Assignseg/conv3/weights/read:02.seg/conv3/weights/Initializer/random_uniform:0
l
seg/conv3/biases:0seg/conv3/biases/Assignseg/conv3/biases/read:02$seg/conv3/biases/Initializer/Const:0
p
seg/conv3/bn/beta:0seg/conv3/bn/beta/Assignseg/conv3/bn/beta/read:02%seg/conv3/bn/beta/Initializer/zeros:0
s
seg/conv3/bn/gamma:0seg/conv3/bn/gamma/Assignseg/conv3/bn/gamma/read:02%seg/conv3/bn/gamma/Initializer/ones:0
y
seg/conv5/weights:0seg/conv5/weights/Assignseg/conv5/weights/read:02.seg/conv5/weights/Initializer/random_uniform:0
l
seg/conv5/biases:0seg/conv5/biases/Assignseg/conv5/biases/read:02$seg/conv5/biases/Initializer/Const:0"Ã
losses¸
µ
-layerfilter0_newfea_conv_head_0/weight_loss:0
$layerfilter0_edgefea_0/weight_loss:0
/layerfilter0_self_att_conv_head_0/weight_loss:0
/layerfilter0_neib_att_conv_head_0/weight_loss:0
gapnet00/weight_loss:0
gapnet01/weight_loss:0
-layerfilter1_newfea_conv_head_0/weight_loss:0
$layerfilter1_edgefea_0/weight_loss:0
/layerfilter1_self_att_conv_head_0/weight_loss:0
/layerfilter1_neib_att_conv_head_0/weight_loss:0
-layerfilter1_newfea_conv_head_1/weight_loss:0
$layerfilter1_edgefea_1/weight_loss:0
/layerfilter1_self_att_conv_head_1/weight_loss:0
/layerfilter1_neib_att_conv_head_1/weight_loss:0
gapnet10/weight_loss:0
gapnet11/weight_loss:0
global_expand/weight_loss:0
agg/weight_loss:0"¹
cond_context§£
ª
1layerfilter0_newfea_conv_head_0/bn/cond/cond_text1layerfilter0_newfea_conv_head_0/bn/cond/pred_id:02layerfilter0_newfea_conv_head_0/bn/cond/switch_t:0 *
Ylayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Xlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Vlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
_layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
alayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Xlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Rlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
[layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Xlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Zlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Xlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
alayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
clayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Zlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Tlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Hlayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
<layerfilter0_newfea_conv_head_0/bn/cond/control_dependency:0
1layerfilter0_newfea_conv_head_0/bn/cond/pred_id:0
2layerfilter0_newfea_conv_head_0/bn/cond/switch_t:0
ulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
playerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
wlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
rlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
4layerfilter0_newfea_conv_head_0/bn/moments/Squeeze:0
6layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1:0Ø
ulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0_layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
4layerfilter0_newfea_conv_head_0/bn/moments/Squeeze:0alayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Í
playerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0Ylayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Ü
wlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0alayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
6layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1:0clayerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1Ñ
rlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0[layerfilter0_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Æ
3layerfilter0_newfea_conv_head_0/bn/cond/cond_text_11layerfilter0_newfea_conv_head_0/bn/cond/pred_id:02layerfilter0_newfea_conv_head_0/bn/cond/switch_f:0*§
>layerfilter0_newfea_conv_head_0/bn/cond/control_dependency_1:0
1layerfilter0_newfea_conv_head_0/bn/cond/pred_id:0
2layerfilter0_newfea_conv_head_0/bn/cond/switch_f:0
Ö
3layerfilter0_newfea_conv_head_0/bn/cond_1/cond_text3layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id:04layerfilter0_newfea_conv_head_0/bn/cond_1/switch_t:0 *±
;layerfilter0_newfea_conv_head_0/bn/cond_1/Identity/Switch:1
4layerfilter0_newfea_conv_head_0/bn/cond_1/Identity:0
=layerfilter0_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1
6layerfilter0_newfea_conv_head_0/bn/cond_1/Identity_1:0
3layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id:0
4layerfilter0_newfea_conv_head_0/bn/cond_1/switch_t:0
4layerfilter0_newfea_conv_head_0/bn/moments/Squeeze:0
6layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1:0s
4layerfilter0_newfea_conv_head_0/bn/moments/Squeeze:0;layerfilter0_newfea_conv_head_0/bn/cond_1/Identity/Switch:1w
6layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1:0=layerfilter0_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1
º
5layerfilter0_newfea_conv_head_0/bn/cond_1/cond_text_13layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id:04layerfilter0_newfea_conv_head_0/bn/cond_1/switch_f:0*
4layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_1:0
4layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_1:1
4layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_2:0
4layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_2:1
3layerfilter0_newfea_conv_head_0/bn/cond_1/pred_id:0
4layerfilter0_newfea_conv_head_0/bn/cond_1/switch_f:0
ulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
wlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0­
ulayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:04layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_1:0¯
wlayerfilter0_newfea_conv_head_0/bn/layerfilter0_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:04layerfilter0_newfea_conv_head_0/bn/cond_1/Switch_2:0
ñ
(layerfilter0_edgefea_0/bn/cond/cond_text(layerfilter0_edgefea_0/bn/cond/pred_id:0)layerfilter0_edgefea_0/bn/cond/switch_t:0 *í
Playerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Olayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Mlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Vlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Xlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Olayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Ilayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Rlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Olayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Qlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Olayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Xlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Zlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Qlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Klayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
?layerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/decay:0
3layerfilter0_edgefea_0/bn/cond/control_dependency:0
(layerfilter0_edgefea_0/bn/cond/pred_id:0
)layerfilter0_edgefea_0/bn/cond/switch_t:0
clayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0
elayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
`layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
+layerfilter0_edgefea_0/bn/moments/Squeeze:0
-layerfilter0_edgefea_0/bn/moments/Squeeze_1:0½
clayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0Vlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
+layerfilter0_edgefea_0/bn/moments/Squeeze:0Xlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1²
^layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0Playerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Á
elayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Xlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
-layerfilter0_edgefea_0/bn/moments/Squeeze_1:0Zlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1¶
`layerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0Rlayerfilter0_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1

*layerfilter0_edgefea_0/bn/cond/cond_text_1(layerfilter0_edgefea_0/bn/cond/pred_id:0)layerfilter0_edgefea_0/bn/cond/switch_f:0*
5layerfilter0_edgefea_0/bn/cond/control_dependency_1:0
(layerfilter0_edgefea_0/bn/cond/pred_id:0
)layerfilter0_edgefea_0/bn/cond/switch_f:0
Ï
*layerfilter0_edgefea_0/bn/cond_1/cond_text*layerfilter0_edgefea_0/bn/cond_1/pred_id:0+layerfilter0_edgefea_0/bn/cond_1/switch_t:0 *Å
2layerfilter0_edgefea_0/bn/cond_1/Identity/Switch:1
+layerfilter0_edgefea_0/bn/cond_1/Identity:0
4layerfilter0_edgefea_0/bn/cond_1/Identity_1/Switch:1
-layerfilter0_edgefea_0/bn/cond_1/Identity_1:0
*layerfilter0_edgefea_0/bn/cond_1/pred_id:0
+layerfilter0_edgefea_0/bn/cond_1/switch_t:0
+layerfilter0_edgefea_0/bn/moments/Squeeze:0
-layerfilter0_edgefea_0/bn/moments/Squeeze_1:0a
+layerfilter0_edgefea_0/bn/moments/Squeeze:02layerfilter0_edgefea_0/bn/cond_1/Identity/Switch:1e
-layerfilter0_edgefea_0/bn/moments/Squeeze_1:04layerfilter0_edgefea_0/bn/cond_1/Identity_1/Switch:1

,layerfilter0_edgefea_0/bn/cond_1/cond_text_1*layerfilter0_edgefea_0/bn/cond_1/pred_id:0+layerfilter0_edgefea_0/bn/cond_1/switch_f:0*
+layerfilter0_edgefea_0/bn/cond_1/Switch_1:0
+layerfilter0_edgefea_0/bn/cond_1/Switch_1:1
+layerfilter0_edgefea_0/bn/cond_1/Switch_2:0
+layerfilter0_edgefea_0/bn/cond_1/Switch_2:1
*layerfilter0_edgefea_0/bn/cond_1/pred_id:0
+layerfilter0_edgefea_0/bn/cond_1/switch_f:0
clayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
elayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
clayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0+layerfilter0_edgefea_0/bn/cond_1/Switch_1:0
elayerfilter0_edgefea_0/bn/layerfilter0_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0+layerfilter0_edgefea_0/bn/cond_1/Switch_2:0

3layerfilter0_self_att_conv_head_0/bn/cond/cond_text3layerfilter0_self_att_conv_head_0/bn/cond/pred_id:04layerfilter0_self_att_conv_head_0/bn/cond/switch_t:0 *ç
[layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter0_self_att_conv_head_0/bn/cond/control_dependency:0
3layerfilter0_self_att_conv_head_0/bn/cond/pred_id:0
4layerfilter0_self_att_conv_head_0/bn/cond/switch_t:0
ylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter0_self_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1:0Þ
ylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
6layerfilter0_self_att_conv_head_0/bn/moments/Squeeze:0clayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ó
tlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1â
{layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¡
8layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1:0elayerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1×
vlayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter0_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Ò
5layerfilter0_self_att_conv_head_0/bn/cond/cond_text_13layerfilter0_self_att_conv_head_0/bn/cond/pred_id:04layerfilter0_self_att_conv_head_0/bn/cond/switch_f:0*­
@layerfilter0_self_att_conv_head_0/bn/cond/control_dependency_1:0
3layerfilter0_self_att_conv_head_0/bn/cond/pred_id:0
4layerfilter0_self_att_conv_head_0/bn/cond/switch_f:0
ô
5layerfilter0_self_att_conv_head_0/bn/cond_1/cond_text5layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id:06layerfilter0_self_att_conv_head_0/bn/cond_1/switch_t:0 *É
=layerfilter0_self_att_conv_head_0/bn/cond_1/Identity/Switch:1
6layerfilter0_self_att_conv_head_0/bn/cond_1/Identity:0
?layerfilter0_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
8layerfilter0_self_att_conv_head_0/bn/cond_1/Identity_1:0
5layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter0_self_att_conv_head_0/bn/cond_1/switch_t:0
6layerfilter0_self_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1:0w
6layerfilter0_self_att_conv_head_0/bn/moments/Squeeze:0=layerfilter0_self_att_conv_head_0/bn/cond_1/Identity/Switch:1{
8layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1:0?layerfilter0_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
à
7layerfilter0_self_att_conv_head_0/bn/cond_1/cond_text_15layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id:06layerfilter0_self_att_conv_head_0/bn/cond_1/switch_f:0*µ
6layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_1:0
6layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_1:1
6layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_2:0
6layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_2:1
5layerfilter0_self_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter0_self_att_conv_head_0/bn/cond_1/switch_f:0
ylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_1:0µ
{layerfilter0_self_att_conv_head_0/bn/layerfilter0_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter0_self_att_conv_head_0/bn/cond_1/Switch_2:0

3layerfilter0_neib_att_conv_head_0/bn/cond/cond_text3layerfilter0_neib_att_conv_head_0/bn/cond/pred_id:04layerfilter0_neib_att_conv_head_0/bn/cond/switch_t:0 *ç
[layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter0_neib_att_conv_head_0/bn/cond/control_dependency:0
3layerfilter0_neib_att_conv_head_0/bn/cond/pred_id:0
4layerfilter0_neib_att_conv_head_0/bn/cond/switch_t:0
ylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1:0Þ
ylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
6layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze:0clayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ó
tlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1â
{layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¡
8layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1:0elayerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1×
vlayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter0_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Ò
5layerfilter0_neib_att_conv_head_0/bn/cond/cond_text_13layerfilter0_neib_att_conv_head_0/bn/cond/pred_id:04layerfilter0_neib_att_conv_head_0/bn/cond/switch_f:0*­
@layerfilter0_neib_att_conv_head_0/bn/cond/control_dependency_1:0
3layerfilter0_neib_att_conv_head_0/bn/cond/pred_id:0
4layerfilter0_neib_att_conv_head_0/bn/cond/switch_f:0
ô
5layerfilter0_neib_att_conv_head_0/bn/cond_1/cond_text5layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id:06layerfilter0_neib_att_conv_head_0/bn/cond_1/switch_t:0 *É
=layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1
6layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity:0
?layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
8layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity_1:0
5layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter0_neib_att_conv_head_0/bn/cond_1/switch_t:0
6layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1:0w
6layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze:0=layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1{
8layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1:0?layerfilter0_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
à
7layerfilter0_neib_att_conv_head_0/bn/cond_1/cond_text_15layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id:06layerfilter0_neib_att_conv_head_0/bn/cond_1/switch_f:0*µ
6layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_1:0
6layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_1:1
6layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_2:0
6layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_2:1
5layerfilter0_neib_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter0_neib_att_conv_head_0/bn/cond_1/switch_f:0
ylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_1:0µ
{layerfilter0_neib_att_conv_head_0/bn/layerfilter0_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter0_neib_att_conv_head_0/bn/cond_1/Switch_2:0
Á
gapnet00/bn/cond/cond_textgapnet00/bn/cond/pred_id:0gapnet00/bn/cond/switch_t:0 *ç
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
gapnet00/bn/moments/Squeeze_1:0
Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:0Hgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1k
gapnet00/bn/moments/Squeeze:0Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Bgapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage:0Bgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Jgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1o
gapnet00/bn/moments/Squeeze_1:0Lgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Dgapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage:0Dgapnet00/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
»
gapnet00/bn/cond/cond_text_1gapnet00/bn/cond/pred_id:0gapnet00/bn/cond/switch_f:0*b
'gapnet00/bn/cond/control_dependency_1:0
gapnet00/bn/cond/pred_id:0
gapnet00/bn/cond/switch_f:0
ý
gapnet00/bn/cond_1/cond_textgapnet00/bn/cond_1/pred_id:0gapnet00/bn/cond_1/switch_t:0 *
$gapnet00/bn/cond_1/Identity/Switch:1
gapnet00/bn/cond_1/Identity:0
&gapnet00/bn/cond_1/Identity_1/Switch:1
gapnet00/bn/cond_1/Identity_1:0
gapnet00/bn/cond_1/pred_id:0
gapnet00/bn/cond_1/switch_t:0
gapnet00/bn/moments/Squeeze:0
gapnet00/bn/moments/Squeeze_1:0E
gapnet00/bn/moments/Squeeze:0$gapnet00/bn/cond_1/Identity/Switch:1I
gapnet00/bn/moments/Squeeze_1:0&gapnet00/bn/cond_1/Identity_1/Switch:1

gapnet00/bn/cond_1/cond_text_1gapnet00/bn/cond_1/pred_id:0gapnet00/bn/cond_1/switch_f:0*£
gapnet00/bn/cond_1/Switch_1:0
gapnet00/bn/cond_1/Switch_1:1
gapnet00/bn/cond_1/Switch_2:0
gapnet00/bn/cond_1/Switch_2:1
gapnet00/bn/cond_1/pred_id:0
gapnet00/bn/cond_1/switch_f:0
Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0h
Ggapnet00/bn/gapnet00/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet00/bn/cond_1/Switch_1:0j
Igapnet00/bn/gapnet00/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet00/bn/cond_1/Switch_2:0
Á
gapnet01/bn/cond/cond_textgapnet01/bn/cond/pred_id:0gapnet01/bn/cond/switch_t:0 *ç
Bgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Agapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
?gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Hgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Jgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Agapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
;gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Dgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Agapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Cgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Agapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Jgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Lgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Cgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
=gapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
1gapnet01/bn/cond/ExponentialMovingAverage/decay:0
%gapnet01/bn/cond/control_dependency:0
gapnet01/bn/cond/pred_id:0
gapnet01/bn/cond/switch_t:0
Ggapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Bgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage:0
Igapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
Dgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage:0
gapnet01/bn/moments/Squeeze:0
gapnet01/bn/moments/Squeeze_1:0
Ggapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/read:0Hgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1k
gapnet01/bn/moments/Squeeze:0Jgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Bgapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage:0Bgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Igapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Jgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1o
gapnet01/bn/moments/Squeeze_1:0Lgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Dgapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage:0Dgapnet01/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
»
gapnet01/bn/cond/cond_text_1gapnet01/bn/cond/pred_id:0gapnet01/bn/cond/switch_f:0*b
'gapnet01/bn/cond/control_dependency_1:0
gapnet01/bn/cond/pred_id:0
gapnet01/bn/cond/switch_f:0
ý
gapnet01/bn/cond_1/cond_textgapnet01/bn/cond_1/pred_id:0gapnet01/bn/cond_1/switch_t:0 *
$gapnet01/bn/cond_1/Identity/Switch:1
gapnet01/bn/cond_1/Identity:0
&gapnet01/bn/cond_1/Identity_1/Switch:1
gapnet01/bn/cond_1/Identity_1:0
gapnet01/bn/cond_1/pred_id:0
gapnet01/bn/cond_1/switch_t:0
gapnet01/bn/moments/Squeeze:0
gapnet01/bn/moments/Squeeze_1:0E
gapnet01/bn/moments/Squeeze:0$gapnet01/bn/cond_1/Identity/Switch:1I
gapnet01/bn/moments/Squeeze_1:0&gapnet01/bn/cond_1/Identity_1/Switch:1

gapnet01/bn/cond_1/cond_text_1gapnet01/bn/cond_1/pred_id:0gapnet01/bn/cond_1/switch_f:0*£
gapnet01/bn/cond_1/Switch_1:0
gapnet01/bn/cond_1/Switch_1:1
gapnet01/bn/cond_1/Switch_2:0
gapnet01/bn/cond_1/Switch_2:1
gapnet01/bn/cond_1/pred_id:0
gapnet01/bn/cond_1/switch_f:0
Ggapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Igapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0h
Ggapnet01/bn/gapnet01/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet01/bn/cond_1/Switch_1:0j
Igapnet01/bn/gapnet01/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet01/bn/cond_1/Switch_2:0
ª
1layerfilter1_newfea_conv_head_0/bn/cond/cond_text1layerfilter1_newfea_conv_head_0/bn/cond/pred_id:02layerfilter1_newfea_conv_head_0/bn/cond/switch_t:0 *
Ylayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Xlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Vlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
_layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
alayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Xlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Rlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
[layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Xlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Zlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Xlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
alayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
clayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Zlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Tlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Hlayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
<layerfilter1_newfea_conv_head_0/bn/cond/control_dependency:0
1layerfilter1_newfea_conv_head_0/bn/cond/pred_id:0
2layerfilter1_newfea_conv_head_0/bn/cond/switch_t:0
ulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
playerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
wlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
rlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
4layerfilter1_newfea_conv_head_0/bn/moments/Squeeze:0
6layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1:0Ø
ulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0_layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
4layerfilter1_newfea_conv_head_0/bn/moments/Squeeze:0alayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Í
playerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0Ylayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Ü
wlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0alayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
6layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1:0clayerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1Ñ
rlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0[layerfilter1_newfea_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Æ
3layerfilter1_newfea_conv_head_0/bn/cond/cond_text_11layerfilter1_newfea_conv_head_0/bn/cond/pred_id:02layerfilter1_newfea_conv_head_0/bn/cond/switch_f:0*§
>layerfilter1_newfea_conv_head_0/bn/cond/control_dependency_1:0
1layerfilter1_newfea_conv_head_0/bn/cond/pred_id:0
2layerfilter1_newfea_conv_head_0/bn/cond/switch_f:0
Ö
3layerfilter1_newfea_conv_head_0/bn/cond_1/cond_text3layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id:04layerfilter1_newfea_conv_head_0/bn/cond_1/switch_t:0 *±
;layerfilter1_newfea_conv_head_0/bn/cond_1/Identity/Switch:1
4layerfilter1_newfea_conv_head_0/bn/cond_1/Identity:0
=layerfilter1_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1
6layerfilter1_newfea_conv_head_0/bn/cond_1/Identity_1:0
3layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id:0
4layerfilter1_newfea_conv_head_0/bn/cond_1/switch_t:0
4layerfilter1_newfea_conv_head_0/bn/moments/Squeeze:0
6layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1:0s
4layerfilter1_newfea_conv_head_0/bn/moments/Squeeze:0;layerfilter1_newfea_conv_head_0/bn/cond_1/Identity/Switch:1w
6layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1:0=layerfilter1_newfea_conv_head_0/bn/cond_1/Identity_1/Switch:1
º
5layerfilter1_newfea_conv_head_0/bn/cond_1/cond_text_13layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id:04layerfilter1_newfea_conv_head_0/bn/cond_1/switch_f:0*
4layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_1:0
4layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_1:1
4layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_2:0
4layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_2:1
3layerfilter1_newfea_conv_head_0/bn/cond_1/pred_id:0
4layerfilter1_newfea_conv_head_0/bn/cond_1/switch_f:0
ulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
wlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0­
ulayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:04layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_1:0¯
wlayerfilter1_newfea_conv_head_0/bn/layerfilter1_newfea_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:04layerfilter1_newfea_conv_head_0/bn/cond_1/Switch_2:0
ñ
(layerfilter1_edgefea_0/bn/cond/cond_text(layerfilter1_edgefea_0/bn/cond/pred_id:0)layerfilter1_edgefea_0/bn/cond/switch_t:0 *í
Playerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Olayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Mlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Vlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Xlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Olayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Ilayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Rlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Olayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Qlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Olayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Xlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Zlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Qlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Klayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
?layerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/decay:0
3layerfilter1_edgefea_0/bn/cond/control_dependency:0
(layerfilter1_edgefea_0/bn/cond/pred_id:0
)layerfilter1_edgefea_0/bn/cond/switch_t:0
clayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0
elayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
`layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
+layerfilter1_edgefea_0/bn/moments/Squeeze:0
-layerfilter1_edgefea_0/bn/moments/Squeeze_1:0½
clayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0Vlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
+layerfilter1_edgefea_0/bn/moments/Squeeze:0Xlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1²
^layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage:0Playerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Á
elayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Xlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
-layerfilter1_edgefea_0/bn/moments/Squeeze_1:0Zlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1¶
`layerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0Rlayerfilter1_edgefea_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1

*layerfilter1_edgefea_0/bn/cond/cond_text_1(layerfilter1_edgefea_0/bn/cond/pred_id:0)layerfilter1_edgefea_0/bn/cond/switch_f:0*
5layerfilter1_edgefea_0/bn/cond/control_dependency_1:0
(layerfilter1_edgefea_0/bn/cond/pred_id:0
)layerfilter1_edgefea_0/bn/cond/switch_f:0
Ï
*layerfilter1_edgefea_0/bn/cond_1/cond_text*layerfilter1_edgefea_0/bn/cond_1/pred_id:0+layerfilter1_edgefea_0/bn/cond_1/switch_t:0 *Å
2layerfilter1_edgefea_0/bn/cond_1/Identity/Switch:1
+layerfilter1_edgefea_0/bn/cond_1/Identity:0
4layerfilter1_edgefea_0/bn/cond_1/Identity_1/Switch:1
-layerfilter1_edgefea_0/bn/cond_1/Identity_1:0
*layerfilter1_edgefea_0/bn/cond_1/pred_id:0
+layerfilter1_edgefea_0/bn/cond_1/switch_t:0
+layerfilter1_edgefea_0/bn/moments/Squeeze:0
-layerfilter1_edgefea_0/bn/moments/Squeeze_1:0a
+layerfilter1_edgefea_0/bn/moments/Squeeze:02layerfilter1_edgefea_0/bn/cond_1/Identity/Switch:1e
-layerfilter1_edgefea_0/bn/moments/Squeeze_1:04layerfilter1_edgefea_0/bn/cond_1/Identity_1/Switch:1

,layerfilter1_edgefea_0/bn/cond_1/cond_text_1*layerfilter1_edgefea_0/bn/cond_1/pred_id:0+layerfilter1_edgefea_0/bn/cond_1/switch_f:0*
+layerfilter1_edgefea_0/bn/cond_1/Switch_1:0
+layerfilter1_edgefea_0/bn/cond_1/Switch_1:1
+layerfilter1_edgefea_0/bn/cond_1/Switch_2:0
+layerfilter1_edgefea_0/bn/cond_1/Switch_2:1
*layerfilter1_edgefea_0/bn/cond_1/pred_id:0
+layerfilter1_edgefea_0/bn/cond_1/switch_f:0
clayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
elayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
clayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0+layerfilter1_edgefea_0/bn/cond_1/Switch_1:0
elayerfilter1_edgefea_0/bn/layerfilter1_edgefea_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0+layerfilter1_edgefea_0/bn/cond_1/Switch_2:0

3layerfilter1_self_att_conv_head_0/bn/cond/cond_text3layerfilter1_self_att_conv_head_0/bn/cond/pred_id:04layerfilter1_self_att_conv_head_0/bn/cond/switch_t:0 *ç
[layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter1_self_att_conv_head_0/bn/cond/control_dependency:0
3layerfilter1_self_att_conv_head_0/bn/cond/pred_id:0
4layerfilter1_self_att_conv_head_0/bn/cond/switch_t:0
ylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter1_self_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1:0Þ
ylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
6layerfilter1_self_att_conv_head_0/bn/moments/Squeeze:0clayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ó
tlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1â
{layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¡
8layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1:0elayerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1×
vlayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter1_self_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Ò
5layerfilter1_self_att_conv_head_0/bn/cond/cond_text_13layerfilter1_self_att_conv_head_0/bn/cond/pred_id:04layerfilter1_self_att_conv_head_0/bn/cond/switch_f:0*­
@layerfilter1_self_att_conv_head_0/bn/cond/control_dependency_1:0
3layerfilter1_self_att_conv_head_0/bn/cond/pred_id:0
4layerfilter1_self_att_conv_head_0/bn/cond/switch_f:0
ô
5layerfilter1_self_att_conv_head_0/bn/cond_1/cond_text5layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id:06layerfilter1_self_att_conv_head_0/bn/cond_1/switch_t:0 *É
=layerfilter1_self_att_conv_head_0/bn/cond_1/Identity/Switch:1
6layerfilter1_self_att_conv_head_0/bn/cond_1/Identity:0
?layerfilter1_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
8layerfilter1_self_att_conv_head_0/bn/cond_1/Identity_1:0
5layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter1_self_att_conv_head_0/bn/cond_1/switch_t:0
6layerfilter1_self_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1:0w
6layerfilter1_self_att_conv_head_0/bn/moments/Squeeze:0=layerfilter1_self_att_conv_head_0/bn/cond_1/Identity/Switch:1{
8layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1:0?layerfilter1_self_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
à
7layerfilter1_self_att_conv_head_0/bn/cond_1/cond_text_15layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id:06layerfilter1_self_att_conv_head_0/bn/cond_1/switch_f:0*µ
6layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_1:0
6layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_1:1
6layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_2:0
6layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_2:1
5layerfilter1_self_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter1_self_att_conv_head_0/bn/cond_1/switch_f:0
ylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_1:0µ
{layerfilter1_self_att_conv_head_0/bn/layerfilter1_self_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter1_self_att_conv_head_0/bn/cond_1/Switch_2:0

3layerfilter1_neib_att_conv_head_0/bn/cond/cond_text3layerfilter1_neib_att_conv_head_0/bn/cond/pred_id:04layerfilter1_neib_att_conv_head_0/bn/cond/switch_t:0 *ç
[layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter1_neib_att_conv_head_0/bn/cond/control_dependency:0
3layerfilter1_neib_att_conv_head_0/bn/cond/pred_id:0
4layerfilter1_neib_att_conv_head_0/bn/cond/switch_t:0
ylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1:0Þ
ylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
6layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze:0clayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ó
tlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1â
{layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¡
8layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1:0elayerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1×
vlayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter1_neib_att_conv_head_0/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Ò
5layerfilter1_neib_att_conv_head_0/bn/cond/cond_text_13layerfilter1_neib_att_conv_head_0/bn/cond/pred_id:04layerfilter1_neib_att_conv_head_0/bn/cond/switch_f:0*­
@layerfilter1_neib_att_conv_head_0/bn/cond/control_dependency_1:0
3layerfilter1_neib_att_conv_head_0/bn/cond/pred_id:0
4layerfilter1_neib_att_conv_head_0/bn/cond/switch_f:0
ô
5layerfilter1_neib_att_conv_head_0/bn/cond_1/cond_text5layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id:06layerfilter1_neib_att_conv_head_0/bn/cond_1/switch_t:0 *É
=layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1
6layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity:0
?layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
8layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity_1:0
5layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter1_neib_att_conv_head_0/bn/cond_1/switch_t:0
6layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze:0
8layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1:0w
6layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze:0=layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity/Switch:1{
8layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1:0?layerfilter1_neib_att_conv_head_0/bn/cond_1/Identity_1/Switch:1
à
7layerfilter1_neib_att_conv_head_0/bn/cond_1/cond_text_15layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id:06layerfilter1_neib_att_conv_head_0/bn/cond_1/switch_f:0*µ
6layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_1:0
6layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_1:1
6layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_2:0
6layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_2:1
5layerfilter1_neib_att_conv_head_0/bn/cond_1/pred_id:0
6layerfilter1_neib_att_conv_head_0/bn/cond_1/switch_f:0
ylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_1:0µ
{layerfilter1_neib_att_conv_head_0/bn/layerfilter1_neib_att_conv_head_0/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter1_neib_att_conv_head_0/bn/cond_1/Switch_2:0
ª
1layerfilter1_newfea_conv_head_1/bn/cond/cond_text1layerfilter1_newfea_conv_head_1/bn/cond/pred_id:02layerfilter1_newfea_conv_head_1/bn/cond/switch_t:0 *
Ylayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Xlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Vlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
_layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
alayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Xlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Rlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
[layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Xlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Zlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Xlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
alayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
clayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Zlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Tlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Hlayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/decay:0
<layerfilter1_newfea_conv_head_1/bn/cond/control_dependency:0
1layerfilter1_newfea_conv_head_1/bn/cond/pred_id:0
2layerfilter1_newfea_conv_head_1/bn/cond/switch_t:0
ulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
playerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0
wlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
rlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0
4layerfilter1_newfea_conv_head_1/bn/moments/Squeeze:0
6layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1:0Ø
ulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0_layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
4layerfilter1_newfea_conv_head_1/bn/moments/Squeeze:0alayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Í
playerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0Ylayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Ü
wlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0alayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
6layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1:0clayerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1Ñ
rlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0[layerfilter1_newfea_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Æ
3layerfilter1_newfea_conv_head_1/bn/cond/cond_text_11layerfilter1_newfea_conv_head_1/bn/cond/pred_id:02layerfilter1_newfea_conv_head_1/bn/cond/switch_f:0*§
>layerfilter1_newfea_conv_head_1/bn/cond/control_dependency_1:0
1layerfilter1_newfea_conv_head_1/bn/cond/pred_id:0
2layerfilter1_newfea_conv_head_1/bn/cond/switch_f:0
Ö
3layerfilter1_newfea_conv_head_1/bn/cond_1/cond_text3layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id:04layerfilter1_newfea_conv_head_1/bn/cond_1/switch_t:0 *±
;layerfilter1_newfea_conv_head_1/bn/cond_1/Identity/Switch:1
4layerfilter1_newfea_conv_head_1/bn/cond_1/Identity:0
=layerfilter1_newfea_conv_head_1/bn/cond_1/Identity_1/Switch:1
6layerfilter1_newfea_conv_head_1/bn/cond_1/Identity_1:0
3layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id:0
4layerfilter1_newfea_conv_head_1/bn/cond_1/switch_t:0
4layerfilter1_newfea_conv_head_1/bn/moments/Squeeze:0
6layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1:0s
4layerfilter1_newfea_conv_head_1/bn/moments/Squeeze:0;layerfilter1_newfea_conv_head_1/bn/cond_1/Identity/Switch:1w
6layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1:0=layerfilter1_newfea_conv_head_1/bn/cond_1/Identity_1/Switch:1
º
5layerfilter1_newfea_conv_head_1/bn/cond_1/cond_text_13layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id:04layerfilter1_newfea_conv_head_1/bn/cond_1/switch_f:0*
4layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_1:0
4layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_1:1
4layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_2:0
4layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_2:1
3layerfilter1_newfea_conv_head_1/bn/cond_1/pred_id:0
4layerfilter1_newfea_conv_head_1/bn/cond_1/switch_f:0
ulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
wlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0­
ulayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:04layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_1:0¯
wlayerfilter1_newfea_conv_head_1/bn/layerfilter1_newfea_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:04layerfilter1_newfea_conv_head_1/bn/cond_1/Switch_2:0
ñ
(layerfilter1_edgefea_1/bn/cond/cond_text(layerfilter1_edgefea_1/bn/cond/pred_id:0)layerfilter1_edgefea_1/bn/cond/switch_t:0 *í
Playerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Olayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Mlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Vlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Xlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Olayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Ilayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Rlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Olayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Qlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Olayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Xlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Zlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Qlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Klayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
?layerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/decay:0
3layerfilter1_edgefea_1/bn/cond/control_dependency:0
(layerfilter1_edgefea_1/bn/cond/pred_id:0
)layerfilter1_edgefea_1/bn/cond/switch_t:0
clayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage:0
elayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
`layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0
+layerfilter1_edgefea_1/bn/moments/Squeeze:0
-layerfilter1_edgefea_1/bn/moments/Squeeze_1:0½
clayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0Vlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
+layerfilter1_edgefea_1/bn/moments/Squeeze:0Xlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1²
^layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage:0Playerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1Á
elayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Xlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
-layerfilter1_edgefea_1/bn/moments/Squeeze_1:0Zlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1¶
`layerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0Rlayerfilter1_edgefea_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1

*layerfilter1_edgefea_1/bn/cond/cond_text_1(layerfilter1_edgefea_1/bn/cond/pred_id:0)layerfilter1_edgefea_1/bn/cond/switch_f:0*
5layerfilter1_edgefea_1/bn/cond/control_dependency_1:0
(layerfilter1_edgefea_1/bn/cond/pred_id:0
)layerfilter1_edgefea_1/bn/cond/switch_f:0
Ï
*layerfilter1_edgefea_1/bn/cond_1/cond_text*layerfilter1_edgefea_1/bn/cond_1/pred_id:0+layerfilter1_edgefea_1/bn/cond_1/switch_t:0 *Å
2layerfilter1_edgefea_1/bn/cond_1/Identity/Switch:1
+layerfilter1_edgefea_1/bn/cond_1/Identity:0
4layerfilter1_edgefea_1/bn/cond_1/Identity_1/Switch:1
-layerfilter1_edgefea_1/bn/cond_1/Identity_1:0
*layerfilter1_edgefea_1/bn/cond_1/pred_id:0
+layerfilter1_edgefea_1/bn/cond_1/switch_t:0
+layerfilter1_edgefea_1/bn/moments/Squeeze:0
-layerfilter1_edgefea_1/bn/moments/Squeeze_1:0a
+layerfilter1_edgefea_1/bn/moments/Squeeze:02layerfilter1_edgefea_1/bn/cond_1/Identity/Switch:1e
-layerfilter1_edgefea_1/bn/moments/Squeeze_1:04layerfilter1_edgefea_1/bn/cond_1/Identity_1/Switch:1

,layerfilter1_edgefea_1/bn/cond_1/cond_text_1*layerfilter1_edgefea_1/bn/cond_1/pred_id:0+layerfilter1_edgefea_1/bn/cond_1/switch_f:0*
+layerfilter1_edgefea_1/bn/cond_1/Switch_1:0
+layerfilter1_edgefea_1/bn/cond_1/Switch_1:1
+layerfilter1_edgefea_1/bn/cond_1/Switch_2:0
+layerfilter1_edgefea_1/bn/cond_1/Switch_2:1
*layerfilter1_edgefea_1/bn/cond_1/pred_id:0
+layerfilter1_edgefea_1/bn/cond_1/switch_f:0
clayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
elayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
clayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0+layerfilter1_edgefea_1/bn/cond_1/Switch_1:0
elayerfilter1_edgefea_1/bn/layerfilter1_edgefea_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0+layerfilter1_edgefea_1/bn/cond_1/Switch_2:0

3layerfilter1_self_att_conv_head_1/bn/cond/cond_text3layerfilter1_self_att_conv_head_1/bn/cond/pred_id:04layerfilter1_self_att_conv_head_1/bn/cond/switch_t:0 *ç
[layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter1_self_att_conv_head_1/bn/cond/control_dependency:0
3layerfilter1_self_att_conv_head_1/bn/cond/pred_id:0
4layerfilter1_self_att_conv_head_1/bn/cond/switch_t:0
ylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter1_self_att_conv_head_1/bn/moments/Squeeze:0
8layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1:0Þ
ylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
6layerfilter1_self_att_conv_head_1/bn/moments/Squeeze:0clayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ó
tlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1â
{layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¡
8layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1:0elayerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1×
vlayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter1_self_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Ò
5layerfilter1_self_att_conv_head_1/bn/cond/cond_text_13layerfilter1_self_att_conv_head_1/bn/cond/pred_id:04layerfilter1_self_att_conv_head_1/bn/cond/switch_f:0*­
@layerfilter1_self_att_conv_head_1/bn/cond/control_dependency_1:0
3layerfilter1_self_att_conv_head_1/bn/cond/pred_id:0
4layerfilter1_self_att_conv_head_1/bn/cond/switch_f:0
ô
5layerfilter1_self_att_conv_head_1/bn/cond_1/cond_text5layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id:06layerfilter1_self_att_conv_head_1/bn/cond_1/switch_t:0 *É
=layerfilter1_self_att_conv_head_1/bn/cond_1/Identity/Switch:1
6layerfilter1_self_att_conv_head_1/bn/cond_1/Identity:0
?layerfilter1_self_att_conv_head_1/bn/cond_1/Identity_1/Switch:1
8layerfilter1_self_att_conv_head_1/bn/cond_1/Identity_1:0
5layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id:0
6layerfilter1_self_att_conv_head_1/bn/cond_1/switch_t:0
6layerfilter1_self_att_conv_head_1/bn/moments/Squeeze:0
8layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1:0w
6layerfilter1_self_att_conv_head_1/bn/moments/Squeeze:0=layerfilter1_self_att_conv_head_1/bn/cond_1/Identity/Switch:1{
8layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1:0?layerfilter1_self_att_conv_head_1/bn/cond_1/Identity_1/Switch:1
à
7layerfilter1_self_att_conv_head_1/bn/cond_1/cond_text_15layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id:06layerfilter1_self_att_conv_head_1/bn/cond_1/switch_f:0*µ
6layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_1:0
6layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_1:1
6layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_2:0
6layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_2:1
5layerfilter1_self_att_conv_head_1/bn/cond_1/pred_id:0
6layerfilter1_self_att_conv_head_1/bn/cond_1/switch_f:0
ylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_1:0µ
{layerfilter1_self_att_conv_head_1/bn/layerfilter1_self_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter1_self_att_conv_head_1/bn/cond_1/Switch_2:0

3layerfilter1_neib_att_conv_head_1/bn/cond/cond_text3layerfilter1_neib_att_conv_head_1/bn/cond/pred_id:04layerfilter1_neib_att_conv_head_1/bn/cond/switch_t:0 *ç
[layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Zlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Xlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
alayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
clayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Zlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
Tlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
]layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Zlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
\layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Zlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
clayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
elayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
\layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Vlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
Jlayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/decay:0
>layerfilter1_neib_att_conv_head_1/bn/cond/control_dependency:0
3layerfilter1_neib_att_conv_head_1/bn/cond/pred_id:0
4layerfilter1_neib_att_conv_head_1/bn/cond/switch_t:0
ylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
tlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0
{layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
vlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0
6layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze:0
8layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1:0Þ
ylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0alayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
6layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze:0clayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1Ó
tlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage:0[layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1â
{layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0clayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1¡
8layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1:0elayerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1×
vlayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage:0]layerfilter1_neib_att_conv_head_1/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Ò
5layerfilter1_neib_att_conv_head_1/bn/cond/cond_text_13layerfilter1_neib_att_conv_head_1/bn/cond/pred_id:04layerfilter1_neib_att_conv_head_1/bn/cond/switch_f:0*­
@layerfilter1_neib_att_conv_head_1/bn/cond/control_dependency_1:0
3layerfilter1_neib_att_conv_head_1/bn/cond/pred_id:0
4layerfilter1_neib_att_conv_head_1/bn/cond/switch_f:0
ô
5layerfilter1_neib_att_conv_head_1/bn/cond_1/cond_text5layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id:06layerfilter1_neib_att_conv_head_1/bn/cond_1/switch_t:0 *É
=layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity/Switch:1
6layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity:0
?layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity_1/Switch:1
8layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity_1:0
5layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id:0
6layerfilter1_neib_att_conv_head_1/bn/cond_1/switch_t:0
6layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze:0
8layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1:0w
6layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze:0=layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity/Switch:1{
8layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1:0?layerfilter1_neib_att_conv_head_1/bn/cond_1/Identity_1/Switch:1
à
7layerfilter1_neib_att_conv_head_1/bn/cond_1/cond_text_15layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id:06layerfilter1_neib_att_conv_head_1/bn/cond_1/switch_f:0*µ
6layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_1:0
6layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_1:1
6layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_2:0
6layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_2:1
5layerfilter1_neib_att_conv_head_1/bn/cond_1/pred_id:0
6layerfilter1_neib_att_conv_head_1/bn/cond_1/switch_f:0
ylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:0
{layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0³
ylayerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze/ExponentialMovingAverage/read:06layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_1:0µ
{layerfilter1_neib_att_conv_head_1/bn/layerfilter1_neib_att_conv_head_1/bn/moments/Squeeze_1/ExponentialMovingAverage/read:06layerfilter1_neib_att_conv_head_1/bn/cond_1/Switch_2:0
Á
gapnet10/bn/cond/cond_textgapnet10/bn/cond/pred_id:0gapnet10/bn/cond/switch_t:0 *ç
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
gapnet10/bn/moments/Squeeze_1:0
Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:0Hgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1k
gapnet10/bn/moments/Squeeze:0Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Bgapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage:0Bgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Jgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1o
gapnet10/bn/moments/Squeeze_1:0Lgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Dgapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage:0Dgapnet10/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
»
gapnet10/bn/cond/cond_text_1gapnet10/bn/cond/pred_id:0gapnet10/bn/cond/switch_f:0*b
'gapnet10/bn/cond/control_dependency_1:0
gapnet10/bn/cond/pred_id:0
gapnet10/bn/cond/switch_f:0
ý
gapnet10/bn/cond_1/cond_textgapnet10/bn/cond_1/pred_id:0gapnet10/bn/cond_1/switch_t:0 *
$gapnet10/bn/cond_1/Identity/Switch:1
gapnet10/bn/cond_1/Identity:0
&gapnet10/bn/cond_1/Identity_1/Switch:1
gapnet10/bn/cond_1/Identity_1:0
gapnet10/bn/cond_1/pred_id:0
gapnet10/bn/cond_1/switch_t:0
gapnet10/bn/moments/Squeeze:0
gapnet10/bn/moments/Squeeze_1:0E
gapnet10/bn/moments/Squeeze:0$gapnet10/bn/cond_1/Identity/Switch:1I
gapnet10/bn/moments/Squeeze_1:0&gapnet10/bn/cond_1/Identity_1/Switch:1

gapnet10/bn/cond_1/cond_text_1gapnet10/bn/cond_1/pred_id:0gapnet10/bn/cond_1/switch_f:0*£
gapnet10/bn/cond_1/Switch_1:0
gapnet10/bn/cond_1/Switch_1:1
gapnet10/bn/cond_1/Switch_2:0
gapnet10/bn/cond_1/Switch_2:1
gapnet10/bn/cond_1/pred_id:0
gapnet10/bn/cond_1/switch_f:0
Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0h
Ggapnet10/bn/gapnet10/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet10/bn/cond_1/Switch_1:0j
Igapnet10/bn/gapnet10/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet10/bn/cond_1/Switch_2:0
Á
gapnet11/bn/cond/cond_textgapnet11/bn/cond/pred_id:0gapnet11/bn/cond/switch_t:0 *ç
Bgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Agapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
?gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Hgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Jgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Agapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
;gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Dgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Agapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Cgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Agapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Jgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Lgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Cgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
=gapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
1gapnet11/bn/cond/ExponentialMovingAverage/decay:0
%gapnet11/bn/cond/control_dependency:0
gapnet11/bn/cond/pred_id:0
gapnet11/bn/cond/switch_t:0
Ggapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Bgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage:0
Igapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
Dgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage:0
gapnet11/bn/moments/Squeeze:0
gapnet11/bn/moments/Squeeze_1:0
Ggapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/read:0Hgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1k
gapnet11/bn/moments/Squeeze:0Jgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Bgapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage:0Bgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Igapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Jgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1o
gapnet11/bn/moments/Squeeze_1:0Lgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Dgapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage:0Dgapnet11/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
»
gapnet11/bn/cond/cond_text_1gapnet11/bn/cond/pred_id:0gapnet11/bn/cond/switch_f:0*b
'gapnet11/bn/cond/control_dependency_1:0
gapnet11/bn/cond/pred_id:0
gapnet11/bn/cond/switch_f:0
ý
gapnet11/bn/cond_1/cond_textgapnet11/bn/cond_1/pred_id:0gapnet11/bn/cond_1/switch_t:0 *
$gapnet11/bn/cond_1/Identity/Switch:1
gapnet11/bn/cond_1/Identity:0
&gapnet11/bn/cond_1/Identity_1/Switch:1
gapnet11/bn/cond_1/Identity_1:0
gapnet11/bn/cond_1/pred_id:0
gapnet11/bn/cond_1/switch_t:0
gapnet11/bn/moments/Squeeze:0
gapnet11/bn/moments/Squeeze_1:0E
gapnet11/bn/moments/Squeeze:0$gapnet11/bn/cond_1/Identity/Switch:1I
gapnet11/bn/moments/Squeeze_1:0&gapnet11/bn/cond_1/Identity_1/Switch:1

gapnet11/bn/cond_1/cond_text_1gapnet11/bn/cond_1/pred_id:0gapnet11/bn/cond_1/switch_f:0*£
gapnet11/bn/cond_1/Switch_1:0
gapnet11/bn/cond_1/Switch_1:1
gapnet11/bn/cond_1/Switch_2:0
gapnet11/bn/cond_1/Switch_2:1
gapnet11/bn/cond_1/pred_id:0
gapnet11/bn/cond_1/switch_f:0
Ggapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Igapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0h
Ggapnet11/bn/gapnet11/bn/moments/Squeeze/ExponentialMovingAverage/read:0gapnet11/bn/cond_1/Switch_1:0j
Igapnet11/bn/gapnet11/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0gapnet11/bn/cond_1/Switch_2:0
¶
global_expand/bn/cond/cond_textglobal_expand/bn/cond/pred_id:0 global_expand/bn/cond/switch_t:0 *Í
Gglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
Fglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
Dglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Mglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Oglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Fglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
@global_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
Iglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Fglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
Hglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
Fglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Oglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Qglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Hglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
Bglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
6global_expand/bn/cond/ExponentialMovingAverage/decay:0
*global_expand/bn/cond/control_dependency:0
global_expand/bn/cond/pred_id:0
 global_expand/bn/cond/switch_t:0
Qglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Lglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage:0
Sglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
Nglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage:0
"global_expand/bn/moments/Squeeze:0
$global_expand/bn/moments/Squeeze_1:0¢
Qglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/read:0Mglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1u
"global_expand/bn/moments/Squeeze:0Oglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
Lglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage:0Gglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1¦
Sglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Oglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1y
$global_expand/bn/moments/Squeeze_1:0Qglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
Nglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage:0Iglobal_expand/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
Ù
!global_expand/bn/cond/cond_text_1global_expand/bn/cond/pred_id:0 global_expand/bn/cond/switch_f:0*q
,global_expand/bn/cond/control_dependency_1:0
global_expand/bn/cond/pred_id:0
 global_expand/bn/cond/switch_f:0
È
!global_expand/bn/cond_1/cond_text!global_expand/bn/cond_1/pred_id:0"global_expand/bn/cond_1/switch_t:0 *Ù
)global_expand/bn/cond_1/Identity/Switch:1
"global_expand/bn/cond_1/Identity:0
+global_expand/bn/cond_1/Identity_1/Switch:1
$global_expand/bn/cond_1/Identity_1:0
!global_expand/bn/cond_1/pred_id:0
"global_expand/bn/cond_1/switch_t:0
"global_expand/bn/moments/Squeeze:0
$global_expand/bn/moments/Squeeze_1:0O
"global_expand/bn/moments/Squeeze:0)global_expand/bn/cond_1/Identity/Switch:1S
$global_expand/bn/moments/Squeeze_1:0+global_expand/bn/cond_1/Identity_1/Switch:1
â
#global_expand/bn/cond_1/cond_text_1!global_expand/bn/cond_1/pred_id:0"global_expand/bn/cond_1/switch_f:0*ó
"global_expand/bn/cond_1/Switch_1:0
"global_expand/bn/cond_1/Switch_1:1
"global_expand/bn/cond_1/Switch_2:0
"global_expand/bn/cond_1/Switch_2:1
!global_expand/bn/cond_1/pred_id:0
"global_expand/bn/cond_1/switch_f:0
Qglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/read:0
Sglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0w
Qglobal_expand/bn/global_expand/bn/moments/Squeeze/ExponentialMovingAverage/read:0"global_expand/bn/cond_1/Switch_1:0y
Sglobal_expand/bn/global_expand/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0"global_expand/bn/cond_1/Switch_2:0
Ê
agg/bn/cond/cond_textagg/bn/cond/pred_id:0agg/bn/cond/switch_t:0 *ÿ
=agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/read:0
8agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage:0
?agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
:agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage:0
=agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/mul:0
<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub/x:0
:agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub:0
Cagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1
Eagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1
<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1:0
6agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg:0
?agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1
<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/mul:0
>agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub/x:0
<agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub:0
Eagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1
Gagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1
>agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1:0
8agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1:0
,agg/bn/cond/ExponentialMovingAverage/decay:0
 agg/bn/cond/control_dependency:0
agg/bn/cond/pred_id:0
agg/bn/cond/switch_t:0
agg/bn/moments/Squeeze:0
agg/bn/moments/Squeeze_1:0
=agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/read:0Cagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch:1a
agg/bn/moments/Squeeze:0Eagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/sub_1/Switch_1:1y
8agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage:0=agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg/Switch:1
?agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0Eagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch:1e
agg/bn/moments/Squeeze_1:0Gagg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/sub_1/Switch_1:1}
:agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage:0?agg/bn/cond/ExponentialMovingAverage/AssignMovingAvg_1/Switch:1

agg/bn/cond/cond_text_1agg/bn/cond/pred_id:0agg/bn/cond/switch_f:0*S
"agg/bn/cond/control_dependency_1:0
agg/bn/cond/pred_id:0
agg/bn/cond/switch_f:0
²
agg/bn/cond_1/cond_textagg/bn/cond_1/pred_id:0agg/bn/cond_1/switch_t:0 *á
agg/bn/cond_1/Identity/Switch:1
agg/bn/cond_1/Identity:0
!agg/bn/cond_1/Identity_1/Switch:1
agg/bn/cond_1/Identity_1:0
agg/bn/cond_1/pred_id:0
agg/bn/cond_1/switch_t:0
agg/bn/moments/Squeeze:0
agg/bn/moments/Squeeze_1:0;
agg/bn/moments/Squeeze:0agg/bn/cond_1/Identity/Switch:1?
agg/bn/moments/Squeeze_1:0!agg/bn/cond_1/Identity_1/Switch:1
¤
agg/bn/cond_1/cond_text_1agg/bn/cond_1/pred_id:0agg/bn/cond_1/switch_f:0*Ó
=agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/read:0
?agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0
agg/bn/cond_1/Switch_1:0
agg/bn/cond_1/Switch_1:1
agg/bn/cond_1/Switch_2:0
agg/bn/cond_1/Switch_2:1
agg/bn/cond_1/pred_id:0
agg/bn/cond_1/switch_f:0Y
=agg/bn/agg/bn/moments/Squeeze/ExponentialMovingAverage/read:0agg/bn/cond_1/Switch_1:0[
?agg/bn/agg/bn/moments/Squeeze_1/ExponentialMovingAverage/read:0agg/bn/cond_1/Switch_2:0
å
seg/conv2/bn/cond/cond_textseg/conv2/bn/cond/pred_id:0seg/conv2/bn/cond/switch_t:0 *
seg/conv2/BiasAdd:0
seg/conv2/bn/beta/read:0
!seg/conv2/bn/cond/Assign/Switch:1
seg/conv2/bn/cond/Assign:0
#seg/conv2/bn/cond/Assign_1/Switch:1
seg/conv2/bn/cond/Assign_1:0
seg/conv2/bn/cond/add:0
seg/conv2/bn/cond/add_1:0
#seg/conv2/bn/cond/batchnorm/Rsqrt:0
#seg/conv2/bn/cond/batchnorm/add/y:0
!seg/conv2/bn/cond/batchnorm/add:0
#seg/conv2/bn/cond/batchnorm/add_1:0
(seg/conv2/bn/cond/batchnorm/mul/Switch:1
!seg/conv2/bn/cond/batchnorm/mul:0
#seg/conv2/bn/cond/batchnorm/mul_1:0
#seg/conv2/bn/cond/batchnorm/mul_2:0
(seg/conv2/bn/cond/batchnorm/sub/Switch:1
!seg/conv2/bn/cond/batchnorm/sub:0
-seg/conv2/bn/cond/moments/SquaredDifference:0
#seg/conv2/bn/cond/moments/Squeeze:0
%seg/conv2/bn/cond/moments/Squeeze_1:0
(seg/conv2/bn/cond/moments/StopGradient:0
'seg/conv2/bn/cond/moments/mean/Switch:1
2seg/conv2/bn/cond/moments/mean/reduction_indices:0
 seg/conv2/bn/cond/moments/mean:0
6seg/conv2/bn/cond/moments/variance/reduction_indices:0
$seg/conv2/bn/cond/moments/variance:0
seg/conv2/bn/cond/mul/Switch:1
seg/conv2/bn/cond/mul/y:0
seg/conv2/bn/cond/mul:0
seg/conv2/bn/cond/mul_1/y:0
seg/conv2/bn/cond/mul_1:0
 seg/conv2/bn/cond/mul_2/Switch:1
seg/conv2/bn/cond/mul_2/y:0
seg/conv2/bn/cond/mul_2:0
seg/conv2/bn/cond/mul_3/y:0
seg/conv2/bn/cond/mul_3:0
seg/conv2/bn/cond/pred_id:0
seg/conv2/bn/cond/switch_t:0
seg/conv2/bn/gamma/read:0
seg/conv2/bn/pop_mean/read:0
seg/conv2/bn/pop_mean:0
seg/conv2/bn/pop_var/read:0
seg/conv2/bn/pop_var:0>
seg/conv2/BiasAdd:0'seg/conv2/bn/cond/moments/mean/Switch:1>
seg/conv2/bn/pop_mean/read:0seg/conv2/bn/cond/mul/Switch:1<
seg/conv2/bn/pop_mean:0!seg/conv2/bn/cond/Assign/Switch:1?
seg/conv2/bn/pop_var/read:0 seg/conv2/bn/cond/mul_2/Switch:1=
seg/conv2/bn/pop_var:0#seg/conv2/bn/cond/Assign_1/Switch:1E
seg/conv2/bn/gamma/read:0(seg/conv2/bn/cond/batchnorm/mul/Switch:1D
seg/conv2/bn/beta/read:0(seg/conv2/bn/cond/batchnorm/sub/Switch:1
	
seg/conv2/bn/cond/cond_text_1seg/conv2/bn/cond/pred_id:0seg/conv2/bn/cond/switch_f:0*Á
seg/conv2/BiasAdd:0
seg/conv2/bn/beta/read:0
%seg/conv2/bn/cond/batchnorm_1/Rsqrt:0
*seg/conv2/bn/cond/batchnorm_1/add/Switch:0
%seg/conv2/bn/cond/batchnorm_1/add/y:0
#seg/conv2/bn/cond/batchnorm_1/add:0
%seg/conv2/bn/cond/batchnorm_1/add_1:0
*seg/conv2/bn/cond/batchnorm_1/mul/Switch:0
#seg/conv2/bn/cond/batchnorm_1/mul:0
,seg/conv2/bn/cond/batchnorm_1/mul_1/Switch:0
%seg/conv2/bn/cond/batchnorm_1/mul_1:0
,seg/conv2/bn/cond/batchnorm_1/mul_2/Switch:0
%seg/conv2/bn/cond/batchnorm_1/mul_2:0
*seg/conv2/bn/cond/batchnorm_1/sub/Switch:0
#seg/conv2/bn/cond/batchnorm_1/sub:0
seg/conv2/bn/cond/pred_id:0
seg/conv2/bn/cond/switch_f:0
seg/conv2/bn/gamma/read:0
seg/conv2/bn/pop_mean/read:0
seg/conv2/bn/pop_var/read:0I
seg/conv2/bn/pop_var/read:0*seg/conv2/bn/cond/batchnorm_1/add/Switch:0G
seg/conv2/bn/gamma/read:0*seg/conv2/bn/cond/batchnorm_1/mul/Switch:0C
seg/conv2/BiasAdd:0,seg/conv2/bn/cond/batchnorm_1/mul_1/Switch:0L
seg/conv2/bn/pop_mean/read:0,seg/conv2/bn/cond/batchnorm_1/mul_2/Switch:0F
seg/conv2/bn/beta/read:0*seg/conv2/bn/cond/batchnorm_1/sub/Switch:0
¥
seg/dp1/cond/cond_textseg/dp1/cond/pred_id:0seg/dp1/cond/switch_t:0 *×
seg/conv2/Relu:0
seg/dp1/cond/dropout/Floor:0
seg/dp1/cond/dropout/Shape:0
seg/dp1/cond/dropout/add:0
!seg/dp1/cond/dropout/div/Switch:1
seg/dp1/cond/dropout/div:0
 seg/dp1/cond/dropout/keep_prob:0
seg/dp1/cond/dropout/mul:0
3seg/dp1/cond/dropout/random_uniform/RandomUniform:0
)seg/dp1/cond/dropout/random_uniform/max:0
)seg/dp1/cond/dropout/random_uniform/min:0
)seg/dp1/cond/dropout/random_uniform/mul:0
)seg/dp1/cond/dropout/random_uniform/sub:0
%seg/dp1/cond/dropout/random_uniform:0
seg/dp1/cond/pred_id:0
seg/dp1/cond/switch_t:05
seg/conv2/Relu:0!seg/dp1/cond/dropout/div/Switch:1
ð
seg/dp1/cond/cond_text_1seg/dp1/cond/pred_id:0seg/dp1/cond/switch_f:0*¢
seg/conv2/Relu:0
seg/dp1/cond/Switch_1:0
seg/dp1/cond/Switch_1:1
seg/dp1/cond/pred_id:0
seg/dp1/cond/switch_f:0+
seg/conv2/Relu:0seg/dp1/cond/Switch_1:0
å
seg/conv3/bn/cond/cond_textseg/conv3/bn/cond/pred_id:0seg/conv3/bn/cond/switch_t:0 *
seg/conv3/BiasAdd:0
seg/conv3/bn/beta/read:0
!seg/conv3/bn/cond/Assign/Switch:1
seg/conv3/bn/cond/Assign:0
#seg/conv3/bn/cond/Assign_1/Switch:1
seg/conv3/bn/cond/Assign_1:0
seg/conv3/bn/cond/add:0
seg/conv3/bn/cond/add_1:0
#seg/conv3/bn/cond/batchnorm/Rsqrt:0
#seg/conv3/bn/cond/batchnorm/add/y:0
!seg/conv3/bn/cond/batchnorm/add:0
#seg/conv3/bn/cond/batchnorm/add_1:0
(seg/conv3/bn/cond/batchnorm/mul/Switch:1
!seg/conv3/bn/cond/batchnorm/mul:0
#seg/conv3/bn/cond/batchnorm/mul_1:0
#seg/conv3/bn/cond/batchnorm/mul_2:0
(seg/conv3/bn/cond/batchnorm/sub/Switch:1
!seg/conv3/bn/cond/batchnorm/sub:0
-seg/conv3/bn/cond/moments/SquaredDifference:0
#seg/conv3/bn/cond/moments/Squeeze:0
%seg/conv3/bn/cond/moments/Squeeze_1:0
(seg/conv3/bn/cond/moments/StopGradient:0
'seg/conv3/bn/cond/moments/mean/Switch:1
2seg/conv3/bn/cond/moments/mean/reduction_indices:0
 seg/conv3/bn/cond/moments/mean:0
6seg/conv3/bn/cond/moments/variance/reduction_indices:0
$seg/conv3/bn/cond/moments/variance:0
seg/conv3/bn/cond/mul/Switch:1
seg/conv3/bn/cond/mul/y:0
seg/conv3/bn/cond/mul:0
seg/conv3/bn/cond/mul_1/y:0
seg/conv3/bn/cond/mul_1:0
 seg/conv3/bn/cond/mul_2/Switch:1
seg/conv3/bn/cond/mul_2/y:0
seg/conv3/bn/cond/mul_2:0
seg/conv3/bn/cond/mul_3/y:0
seg/conv3/bn/cond/mul_3:0
seg/conv3/bn/cond/pred_id:0
seg/conv3/bn/cond/switch_t:0
seg/conv3/bn/gamma/read:0
seg/conv3/bn/pop_mean/read:0
seg/conv3/bn/pop_mean:0
seg/conv3/bn/pop_var/read:0
seg/conv3/bn/pop_var:0>
seg/conv3/BiasAdd:0'seg/conv3/bn/cond/moments/mean/Switch:1>
seg/conv3/bn/pop_mean/read:0seg/conv3/bn/cond/mul/Switch:1<
seg/conv3/bn/pop_mean:0!seg/conv3/bn/cond/Assign/Switch:1?
seg/conv3/bn/pop_var/read:0 seg/conv3/bn/cond/mul_2/Switch:1=
seg/conv3/bn/pop_var:0#seg/conv3/bn/cond/Assign_1/Switch:1E
seg/conv3/bn/gamma/read:0(seg/conv3/bn/cond/batchnorm/mul/Switch:1D
seg/conv3/bn/beta/read:0(seg/conv3/bn/cond/batchnorm/sub/Switch:1
	
seg/conv3/bn/cond/cond_text_1seg/conv3/bn/cond/pred_id:0seg/conv3/bn/cond/switch_f:0*Á
seg/conv3/BiasAdd:0
seg/conv3/bn/beta/read:0
%seg/conv3/bn/cond/batchnorm_1/Rsqrt:0
*seg/conv3/bn/cond/batchnorm_1/add/Switch:0
%seg/conv3/bn/cond/batchnorm_1/add/y:0
#seg/conv3/bn/cond/batchnorm_1/add:0
%seg/conv3/bn/cond/batchnorm_1/add_1:0
*seg/conv3/bn/cond/batchnorm_1/mul/Switch:0
#seg/conv3/bn/cond/batchnorm_1/mul:0
,seg/conv3/bn/cond/batchnorm_1/mul_1/Switch:0
%seg/conv3/bn/cond/batchnorm_1/mul_1:0
,seg/conv3/bn/cond/batchnorm_1/mul_2/Switch:0
%seg/conv3/bn/cond/batchnorm_1/mul_2:0
*seg/conv3/bn/cond/batchnorm_1/sub/Switch:0
#seg/conv3/bn/cond/batchnorm_1/sub:0
seg/conv3/bn/cond/pred_id:0
seg/conv3/bn/cond/switch_f:0
seg/conv3/bn/gamma/read:0
seg/conv3/bn/pop_mean/read:0
seg/conv3/bn/pop_var/read:0I
seg/conv3/bn/pop_var/read:0*seg/conv3/bn/cond/batchnorm_1/add/Switch:0G
seg/conv3/bn/gamma/read:0*seg/conv3/bn/cond/batchnorm_1/mul/Switch:0C
seg/conv3/BiasAdd:0,seg/conv3/bn/cond/batchnorm_1/mul_1/Switch:0L
seg/conv3/bn/pop_mean/read:0,seg/conv3/bn/cond/batchnorm_1/mul_2/Switch:0F
seg/conv3/bn/beta/read:0*seg/conv3/bn/cond/batchnorm_1/sub/Switch:0
¥
seg/dp2/cond/cond_textseg/dp2/cond/pred_id:0seg/dp2/cond/switch_t:0 *×
seg/conv3/Relu:0
seg/dp2/cond/dropout/Floor:0
seg/dp2/cond/dropout/Shape:0
seg/dp2/cond/dropout/add:0
!seg/dp2/cond/dropout/div/Switch:1
seg/dp2/cond/dropout/div:0
 seg/dp2/cond/dropout/keep_prob:0
seg/dp2/cond/dropout/mul:0
3seg/dp2/cond/dropout/random_uniform/RandomUniform:0
)seg/dp2/cond/dropout/random_uniform/max:0
)seg/dp2/cond/dropout/random_uniform/min:0
)seg/dp2/cond/dropout/random_uniform/mul:0
)seg/dp2/cond/dropout/random_uniform/sub:0
%seg/dp2/cond/dropout/random_uniform:0
seg/dp2/cond/pred_id:0
seg/dp2/cond/switch_t:05
seg/conv3/Relu:0!seg/dp2/cond/dropout/div/Switch:1
ð
seg/dp2/cond/cond_text_1seg/dp2/cond/pred_id:0seg/dp2/cond/switch_f:0*¢
seg/conv3/Relu:0
seg/dp2/cond/Switch_1:0
seg/dp2/cond/Switch_1:1
seg/dp2/cond/pred_id:0
seg/dp2/cond/switch_f:0+
seg/conv3/Relu:0seg/dp2/cond/Switch_1:0
¶
cond/cond_textcond/pred_id:0cond/switch_t:0 *
cond/Switch_1:0
cond/Switch_1:1
cond/pred_id:0
cond/switch_t:0
seg/conv5/BiasAdd:0&
seg/conv5/BiasAdd:0cond/Switch_1:1

cond/cond_text_1cond/pred_id:0cond/switch_f:0*â
cond/Rank:0
cond/Reshape/Switch:0
cond/Reshape:0
cond/Reshape_1:0
cond/Shape:0
cond/Shape_1:0
cond/Slice/begin:0
cond/Slice/size:0
cond/Slice:0
cond/Softmax:0
cond/Sub/y:0

cond/Sub:0
cond/concat/axis:0
cond/concat/values_0:0
cond/concat:0
cond/pred_id:0
cond/switch_f:0
seg/conv5/BiasAdd:0,
seg/conv5/BiasAdd:0cond/Reshape/Switch:0"Ù
model_variablesÅÂ
d
BiasAdd/biases:0BiasAdd/biases/AssignBiasAdd/biases/read:02"BiasAdd/biases/Initializer/zeros:0
l
BiasAdd_1/biases:0BiasAdd_1/biases/AssignBiasAdd_1/biases/read:02$BiasAdd_1/biases/Initializer/zeros:0
l
BiasAdd_2/biases:0BiasAdd_2/biases/AssignBiasAdd_2/biases/read:02$BiasAdd_2/biases/Initializer/zeros:0