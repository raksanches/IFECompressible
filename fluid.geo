p0 = newp; Point(p0) = {0.000000, 0.000000, 0.0, 0.010000}; Physical Point('p0') = {p0};
//
p1 = newp; Point(p1) = {0.600000, 0.000000, 0.0, 0.010000}; Physical Point('p1') = {p1};
//
p2 = newp; Point(p2) = {0.600000, 0.200000, 0.0, 0.010000}; Physical Point('p2') = {p2};
//
p3 = newp; Point(p3) = {3.000000, 0.200000, 0.0, 0.010000}; Physical Point('p3') = {p3};
//
p4 = newp; Point(p4) = {3.000000, 1.000000, 0.0, 0.010000}; Physical Point('p4') = {p4};
//
p5 = newp; Point(p5) = {0.000000, 1.000000, 0.0, 0.010000}; Physical Point('p5') = {p5};
//
l0 = newl; Line(l0) = {p0, p1}; Physical Line('l0') = {l0};
//
l1 = newl; Line(l1) = {p1, p2}; Physical Line('l1') = {l1};
//
l2 = newl; Line(l2) = {p2, p3}; Physical Line('l2') = {l2};
//
l3 = newl; Line(l3) = {p3, p4}; Physical Line('l3') = {l3};
//
l4 = newl; Line(l4) = {p4, p5}; Physical Line('l4') = {l4};
//
l5 = newl; Line(l5) = {p5, p0}; Physical Line('l5') = {l5};
//
ll0 = newll; Line Loop(ll0) = {l0, l1, l2, l3, l4, l5};
//
s0 = news; Plane Surface(s0) = {ll0}; Physical Surface('s0') = {s0};
//
Transfinite Line {l0} = 12 Using Progression 1;
//
Transfinite Line {l1} = 6 Using Progression 1;
//
Transfinite Line {l2} = 48 Using Progression 1;
//
Transfinite Line {l3} = 16 Using Progression 1;
//
Transfinite Line {l4} = 60 Using Progression 1;
//
Transfinite Line {l5} = 20 Using Progression 1;
//
