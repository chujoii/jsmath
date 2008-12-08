/* Javascript approximation library for jQuery, v. 0.5.
 *
 * Released under the MIT license by chujoii, December 2008.
 *
 */

function GramMatrix(data, degreeOFApolynomial)
{
	/*
	  Gram matrix for polinom (a + bx + bx^2 + cx^3 + ... zx^n = y)  :

	  sum(xi^0) + sum(xi^1)     + sum(xi^2)     + ... + sum(xi^n)     = sum(yi*xi^0)
	  sum(xi^1) + sum(xi^2)     + sum(xi^3)     + ... + sum(xi^(n+1)) = sum(yi*xi^1)
	  sum(xi^2) + sum(xi^3)     + sum(xi^4)     + ... + sum(xi^(n+2)) = sum(yi*xi^2)
	  ...
	  sum(xi^n) + sum(xi^(n+1)) + sum(xi^(n+2)) + ... + sum(xi^(n+n)) = sum(yi*xi^n)



	  many element repeated:
	  a b c d  y1
	  b c d e  y2
	  c d e f  y3

	*/
	

	var i,j;
	var m = degreeOFApolynomial+1;   // a + bx + cx^2    wc(a,b,c) = polinumLevel+1
	var n = data.length;
	
	var q = new Array(); // for power of points
	
	var a = new Array(); // result
	
	for (i = 0; i < n; ++i)
		q[i] = 1.0;

	for (i = 0; i < m; ++i){
		a[i] = new Array();
		for (j = 0; j < m; ++j)
			a[i].push(0.0);
	}
	a[0][0] = 1.0;
	


	
	var s; // accumulates values for the top row
	var r; // accumulates values for the last column
		
	// generate top row and    rigth column sum(y*x^n)
	for (i = 0; i < m; ++i) {
		s = 0.0;
		r = 0.0;
		for (j = 0; j < n; ++j) {
			s += q[j];             // s = sum(q)
			r += q[j]*data[j][1];  // r = sum(y*q)
			q[j] *= data[j][0]; // q = x^n
		}
		a[0][i] = s;
		a[i][m] = r;
	}
		
	// generate other element of Gram matrix 
	for (i = 1; i < m; ++i) {
		s = 0;		
		for (j = 0; j < n; ++j) {
			s = s + q[j];
			q[j] = q[j]*data[j][0]; // generate   (right-1) column    q = x^n
			if (j<m-1)                     // very strange @%@
				a[i][j] = a[i-1][j+1]; // repeat element
		}
		a[i][m-1] = s;
	}

	return a;
}



function Gaussian_elimination(matrix)
{

	var n = matrix.length;
	var k;
	var i, j;
	var tmp;
	var solution = [];

	//for (k = 0; k < n; ++k)
	//	solution[k] = 0.0;

	//            first part (Forward Elimination)
	for (k = 0; k < n; ++k) {
		tmp = matrix[k][k];
		for (j = k; j < n+1; ++j)
			matrix[k][j] = matrix[k][j] / tmp;
		for (i = k+1; i < n; ++i){
			tmp = matrix[i][k];
			for (j = k; j < n+1; ++j) {
				matrix[i][j] = matrix[i][j] - matrix[k][j]*tmp;
			}
		}
	}

	//            second step back substitution
	solution[n-1] = matrix[n-1][n]/ matrix[n-1][n-1];
	for (k = n-2; k >= 0; --k) {
		tmp  =  matrix[k][n];
		for (j = k+1; j < n; ++j)
			tmp = tmp - solution[j]*matrix[k][j];
		solution[k] = tmp;
	}
	



	return solution;
}



function GramMatrixForParabola(data)
{
	var 
		a11, a12, a13, b1,
		a21, a22, a23, b2,
		a31, a32, a33, b3;
	var x, y;
	var i;

	a11 = 0.0;
	a12 = 0.0;
	a13 = 0.0;
	b1  = 0.0;
	
	a21 = 0.0;
	a22 = 0.0;
	a23 = 0.0;
	b2  = 0.0;
	
	a31 = 0.0;
	a32 = 0.0;
	a33 = 0.0;
	b3  = 0.0;

	for (i = 0; i < data.length; ++i) {
		x = data[i][0];
		y = data[i][1];
		a11 += 1.0;      a12 += x;        a13 += x*x;     b1 += y;
		/* a21 += x;      a22 += x*x;   */ a23 += x*x*x;   b2 += y*x;
		/* a31 += x*x;    a32 += x*x*x; */ a33 += x*x*x*x; b3 += y*x*x;
	}
	a21 = a12;
	a22 = a13;
	a31 = a13;
	a32 = a23;
	
	// matrix solution:
	var d, d1, d2, d3;

	d = 0.0;
	d1 = 0.0;
	d2 = 0.0;
	d3 = 0.0;

	d = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;

	d1 = b1*a22*a33 + a12*a23*b3  + a13*b2*a32  - b1*a23*a32  - a12*b2*a33  - a13*a22*b3;
	d2 = a11*b2*a33 + b1*a23*a31  + a13*a21*b3  - a11*a23*b3  - b1*a21*a33  - a13*b2*a31;
	d3 = a11*a22*b3 + a12*b2*a31  + b1*a21*a32  - a11*b2*a32  - a12*a21*b3  - b1*a22*a31;

	return [d1/d, d2/d, d3/d];
}



function paraboloidpeak(data, searchx, peakwidth)
{

	// use only if (degreeOFApolynomial == 2):    y = a + bx + cx^2    
	var area = [];
	var cur = [];
	var i = 0;
	do{
		cur = data[i];
		if (cur != null){
			if (cur[0] >= searchx-peakwidth)
				area.push(cur);
		}
		else {
			cur = [];
			cur[0] = searchx+peakwidth - 1;		
		}
		i++;
	} while ((cur[0] < searchx+peakwidth) && (i<data.length));

	
	var a, b, c;
	var matrix = [];
	matrix = Gaussian_elimination(GramMatrix(area, 2));
	//matrix = GramMatrixForParabola(area);
	var x = -b/(2.0*c);
	var y = a + b*x + c*x*x;

	return [x, y, a, b, c];
}



function apxFunc (coeff, xmin, xmax, num){
	var x, y;
	var powmax = coeff.length;
	var p;
	var data = new Array();
	var xmaxlimit = xmax + (xmax-xmin)/num; // add last point
	for (var i = xmin; i < xmaxlimit; i += (xmax-xmin)/num) {
		y = 0.0;
		x = 1.0;
		for (p = 0; p < powmax; p += 1){
			y += coeff[p]*x;
			x *=i;
		}
		data.push([i, y]);
	}
	return data;
}


function approximation(data, degreeOFApolynomial, pointsnum)
{
	var a, b, c;
	var xmin, xmax;
	
	var start = data[0];
	var end = data[data.length-1];
	return apxFunc (Gaussian_elimination(GramMatrix(data, degreeOFApolynomial)), start[0], end[0], pointsnum); 
}



function autoYmean(data)
{
	var sum = 0.0;
	var cur;
	for (var i = 0; i < data.length; ++i) {
		cur = data[i];
		if (cur != null)
			sum += cur[1];
	}
	return sum / data.length;
}



function simpleSeachPeak(data, yNoise)
{
	var z;
	var openArea = false;
	var y;
	var listMax = new Array();
	var MaxX;
	var MaxY;
	for (var i = 0; i < data.length; ++i) {
		cur = data[i];
		if (cur != null){
			x = cur[0];
			y = cur[1];
			if (openArea && (y >= yNoise))
				if (MaxY<y) {MaxX = x; MaxY = y;}
			if (!openArea && (y >= yNoise)) {
				openArea = true;
				MaxX = x; MaxY = y;
			}
			if  (openArea && (y < yNoise)){
				openArea = false;
				listMax.push([MaxX, MaxY]);
			}
		}
	}
	return listMax;
}



function searchPeak(data, yNoise, peakwidth)
{
	var listMax = simpleSeachPeak(data, yNoise);
	var cur = [];
	var peakx, peaky, a, b, c;
	var peakcoord = [];
	var listPeaks = new Array();
	for (var i = 0; i < listMax.length; ++i) {
		cur = listMax[i];
		peakcoord = paraboloidpeak(data, cur[0], peakwidth);
		listPeaks.push([peakcoord[0], peakcoord[1]]);
	}
	return listPeaks;
}



