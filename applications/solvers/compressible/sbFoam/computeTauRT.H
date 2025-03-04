volScalarField ux = U.component(vector::X);
volScalarField uy = U.component(vector::Y);

volTensorField gradU = fvc::grad(U);
volScalarField dux_dx = gradU.component(tensor::XX);
volScalarField duy_dx = gradU.component(tensor::XY);
volScalarField dux_dy = gradU.component(tensor::YX);
volScalarField duy_dy = gradU.component(tensor::YY);

volScalarField xCell = mesh.C().component(vector::X);
volScalarField yCell = mesh.C().component(vector::Y);
volScalarField rSqCell = sqr(xCell) + sqr(yCell);

tauRT = muSB/rSqCell * 
		(  sqr(xCell) * duy_dx 
		 + xCell*yCell * ( duy_dy - dux_dx )
		 - sqr(yCell) * dux_dy
		 - xCell * uy
		 + yCell * ux );
