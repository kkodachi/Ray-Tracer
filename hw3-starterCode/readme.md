This program implements all core requirements which are as follows: triangle
intersection, sphere intersection, triangle Phong shading, sphere Phong shading,
shadows rays, and still images. In addition it implements antialiasing and soft
shadows for extra credit. Antialiasing is accomplished by subdividing pixels into
nine sections. Soft shadows are implemented by generating 10 random points inside a
sphere with 0.0-1.0 radius for each light source. The severity can be changed by
multiplying radius by a float but it is currently just x1.0 (line 238). Compile by
calling make in the hw3-starterCode directory. Run by calling "./hw3 <scene of
choice>.scene <outputname>.jpeg". ANTIALIASING and SOFTSHADOWS are booleans that can
turn the effects on and off (line 58). The still images can be found in folders with
the following contents: baseJPEGSS with images both false, softshadowJPEGS with
SOFTSHADOWS true, antialiasingJPEGS with ANTIALIASING true, and combinedJPEGS both
true.