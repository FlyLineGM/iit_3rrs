(* ::Package:: *)

 BeginPackage["rotationUtilsRobotics`"]

rotationUtilsRobotics::usage = "\!\(\*
StyleBox[\"rotationUtilsRobotics\",\nFontWeight->\"Bold\"]\) package contains essential rotational utility modules for robotics.
The package includes the modules \!\(\*
StyleBox[\"DH2T\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"rotx\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"roty\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"rotz\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"mf\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"ms\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"multiplyT\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"DHtable2T\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"axis2skew\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"skew2axis\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"R2J\[Omega]s\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\) and \!\(\*
StyleBox[\"T2Jv\",\nFontWeight->\"Bold\",\nFontSlant->\"Italic\"]\). Use ?\!\(\*
StyleBox[\"ModuleName\",\nFontSlant->\"Italic\"]\) for information on a particular module."

(**************************************************************)
(* Routines to be exported *)
(**************************************************************)
rotx::usage ="\!\(\*
StyleBox[\"rotx\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)] returns the rotation matrix corresponding to rotation by an angle \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\) in CCW sense about x-axis.";

roty::usage ="\!\(\*
StyleBox[\"roty\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)] returns the rotation matrix corresponding to rotation by an angle \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\) in CCW sense about y-axis.";

rotz::usage ="\!\(\*
StyleBox[\"rotz\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)] returns the rotation matrix corresponding to rotation by an angle \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\) in CCW sense about z-axis.";

mf::usage ="\!\(\*
StyleBox[\"mf\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\)] shorthand for MatrixForm of a matrix";

ms::usage = "\!\(\*
StyleBox[\"ms\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\)] shorthand for simplification of a matrix";

DH2T::usage = "\!\(\*
StyleBox[\"DH2T\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[SubscriptBox[\"\[Alpha]\", 
RowBox[{\"i\", \"-\", \"1\"}]],\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[SubscriptBox[\"a\", 
RowBox[{\"i\", \"-\", \"1\"}]],\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[SubscriptBox[\"d\", \"i\"],\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[SubscriptBox[\"\[Theta]\", \"i\"],\nFontSlant->\"Italic\"]\)] returns the transformation matrix from the D-H parameters.";

multiplyT::usage = "\!\(\*
StyleBox[\"multiplyT\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"T1\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"T2\",\nFontSlant->\"Italic\"]\)] multiplies two transformation matrices of dimension 4x4 and returns the resulting matrix.";

DHtable2T::usage = "\!\(\*
StyleBox[\"DHtable2T\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"DHtable\",\nFontSlant->\"Italic\"]\)] returns the 4x4 transformation matrix for a n-link robot, where the input DHtable is the corresponding nx4 matrix.";

axis2skew::usage = "\!\(\*
StyleBox[\"axis2skew\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\)] returns the 3x3 skew-symmetric matrix corresponding to the 3-vector \!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\).";

skew2axis::usage = "\!\(\*
StyleBox[\"skew2axis\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\)] returns the 3-vector corresponding to the 3x3 skew-symmetric \!\(\*
StyleBox[\"matrix\",\nFontSlant->\"Italic\"]\).";

R2J\[Omega]s::usage = "\!\(\*
StyleBox[\"R2J\[Omega]s\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"R\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)] returns the Jacobian corrsponding to the space-fixed angular velocity for a rotation matrix \!\(\*
StyleBox[\"R\",\nFontSlant->\"Italic\"]\) and actuator variable vector \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\).";

T2Jv::usage ="\!\(\*
StyleBox[\"T2Jv\",\nFontWeight->\"Bold\"]\)[\!\(\*
StyleBox[\"T\",\nFontSlant->\"Italic\"]\), \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)] returns the Jacobian which maps the joint rates to the linear velocity."; 

(**************************************************************)
(* Define the routines *)
(**************************************************************)

Begin["`Private`"];


(* rotx[\[Theta]] gives the rotation matrix for CCW rotation by \[Theta] about the X axis *)

rotx[\[Theta]_]:=

	Module[{c\[Theta], s\[Theta]}, (* c\[Theta], s\[Theta] are the local variables used in this module *)

	c\[Theta] = Cos[\[Theta]]; 

	s\[Theta] = Sin[\[Theta]]; 

	Return[{{1,0,0},{0,c\[Theta],-s\[Theta]},{0,s\[Theta],c\[Theta]}}];

]; 


(* roty[\[Theta]] gives the rotation matrix for CCW rotation by \[Theta] about the Y axis *)

roty[\[Theta]_]:=

	Module[{c\[Theta], s\[Theta]}, (* c\[Theta], s\[Theta] are the local variables used in this module *)

	c\[Theta] = Cos[\[Theta]]; 

	s\[Theta] = Sin[\[Theta]]; 

	Return[{{c\[Theta],0,s\[Theta]},{0,1,0},{-s\[Theta],0,c\[Theta]}}];

]; 



(* rotz[\[Theta]] gives the rotation matrix for CCW rotation by \[Theta] about the Z axis *)

rotz[\[Theta]_]:=

	Module[{c\[Theta], s\[Theta]}, (* c\[Theta], s\[Theta] are the local variables used in this module *)

	c\[Theta] = Cos[\[Theta]]; 

	s\[Theta] = Sin[\[Theta]]; 

	Return[{{c\[Theta],-s\[Theta],0},{s\[Theta],c\[Theta],0},{0,0,1}}];

]; 


mf[A_]:= MatrixForm[A]; 

ms[A_]:= MatrixForm[Simplify[A]];



(*DH2T[Subscript[\[Alpha], i-1],Subscript[a, i-1],Subscript[d, i],Subscript[\[Theta], i]] returns the transformation matrix from the D-H parameters.*)

DH2T[{\[Alpha]_,a_,d_,\[Theta]_}]:= Module[{s\[Theta],c\[Theta],s\[Alpha],c\[Alpha],T},
	s\[Theta]=Sin[\[Theta]];
	c\[Theta]=Cos[\[Theta]];
	s\[Alpha]=Sin[\[Alpha]];
	c\[Alpha]=Cos[\[Alpha]];
	T={{c\[Theta],-s\[Theta],0,a},
		{s\[Theta]*c\[Alpha],c\[Theta]*c\[Alpha],-s\[Alpha],-d*s\[Alpha]},
		{s\[Theta]*s\[Alpha],c\[Theta]*s\[Alpha],c\[Alpha],d*c\[Alpha]},
		{0,0,0,1}};
	Return[T];
];




(* multiply two transformation matrices. Instead of directly using the standard built-in matrix multiplication operator, 
we use properties of transformation matrix to develop a computationally efficient module. 
The transformation matrix can be broken down into the rotation matrix and translation vector. 
Rotation matrices are multiplied using standard built-in operator and the displacement vector is multiplied by transpose of R. *)
multiplyT[T1_,T2_] := Module[{R1,R2,d1,d2,T3},
	R1=T1[[1;;3,1;;3]];
	R2=T2[[1;;3,1;;3]];
	d1=T1[[1;;3,4]];
	d2=T2[[1;;3,4]];
	T3=T1;

	T3[[1;;3,1;;3]]=R1.R2;
	T3[[1;;3,4]]=R1.d2+d1;
	Return[T3];
];


(* DHtable2T [DHtable] gives the 4x4 transformation matrix for a n-link robot, where the input DHtable is the corresponding nx4 matrix.*)
DHtable2T[DHtable_]:= Module[{dimrow,dimcol,ii,T},

	DHtable2T::wronginput = " `1` is not a valid DH parameter table! Please inpul a nx4 matrix.";
	(*Validate the input: check if input matrix is nx4*)
	If[Not[MatrixQ[DHtable]]||({dimrow,dimcol}=Dimensions[DHtable];dimcol!=4),
		Message[DHtable2T::wronginput,DHtable];
		Return[];
	];
	(* Iterate till the last frame*)
	T=IdentityMatrix[4];
	For[ii=1,ii<=dimrow,
		T=multiplyT[T,DH2T[DHtable[[ii]]]];
		ii++;
	];
	Return[T];
];




(* axis2skew[d] gives the 3x3 skew-symmetric matrix corresponding to the 3-vector d *)
axis2skew[d_] := Module[{},

	axis2skew::wronginput ="`1` is not a vector of length 3!";

	(*Check to make sure the dimensions are okay*)
	If[!VectorQ[d]||Length[d]!=3,
		Message[axis2skew::wronginput,d];
		Return[];
	];
	(*Return the appropriate matrix*)
	Return[{{0,-d[[3]],d[[2]]},{d[[3]],0,-d[[1]]},{-d[[2]],d[[1]],0}}];
];


(* skew2axis[D] gives the 3-vector corresponding to the 3x3 skew-symmetric matrix D *)
skew2axis[Dm_]:= Module[{},
	skew2axis::wronginput = "`1` is not a 3x3 matrix!";
	(*Check to make sure the dimensions are correct*)
	If[!MatrixQ[Dm]||Dimensions[Dm] !=  {3,3},
	Message[skew2axis::wronginput,Dm];
	Return[];
	]; 

	(* TODO: check for skew-symmetry -- optional *)
	(*Return the appropriate matrix*)
	Return[{Dm[[3,2]],Dm[[1,3]],Dm[[2,1]]}];
];



(* R2J\[Omega]s[R, \[Theta]] gives the Jacobian corrsponding to the space-fixed angular velocity for a rotation R and actuator variable \[Theta]. *)
R2J\[Omega]s[R_, \[Theta]_] := Module[{J\[Omega]s,n, ii}, 
	R2J\[Omega]s::wronginput = "`1` is not a 3x3 matrix!";
	(*Check to make sure the dimensions are correct *)
	If[!MatrixQ[R]||Dimensions[R] != {3,3},
		Message[R2J\[Omega]s::wronginput,R];
		Return[]
	(* TODO: check for orthogonality -- optional *)
	];
	J\[Omega]s = skew2axis[Simplify[D[R, #]. Transpose[R]]]& /@ \[Theta]; 
	(* This means, "skew2axis[Simplify[D[R, #]. Transpose[R]]]&" is a function, with an argument "#", which is then run over the vector \[Theta] by the part "/@ \[Theta]". 
	This is known as "functional programming", and as one can see, it is very compact and powerful.
	Of course, one  can choose to use the traditional "for" loop to do the same thing. *)
	Return[Transpose[J\[Omega]s]];
];



(*T2Jv[T, \[Theta]] returns the Jacobian which maps the joint rates to the linear velocity.*)
T2Jv[T_,\[Theta]_]:=Module[{Jvs,n,ii,d},
	(*Check to make sure the dimensions are correct*)
T2Jv::wronginput = "`1` is not a 4x4 matrix!";
	If[!MatrixQ[T]||Dimensions[T]!={4,4},
		Message[T2Jv::wronginput,T];
		Return[];
	];
	d=T[[1;;3,4]];
	Jvs=Transpose[D[d,#]&/@\[Theta]];
	Return[Jvs];
];


End[];

SetAttributes[rotx,{ReadProtected,Locked}];
SetAttributes[roty,{ReadProtected,Locked}];
SetAttributes[rotz,{ReadProtected,Locked}];
SetAttributes[mf,{ReadProtected,Locked}];
SetAttributes[ms,{ReadProtected,Locked}];
SetAttributes[DHtable2T,{ReadProtected,Locked}];
SetAttributes[T2Jv,{ReadProtected,Locked}];
SetAttributes[R2Jws,{ReadProtected,Locked}];
SetAttributes[skew2axis,{ReadProtected,Locked}];
SetAttributes[axis2skew,{ReadProtected,Locked}];
SetAttributes[multiplyT,{ReadProtected,Locked}];
EndPackage[];
