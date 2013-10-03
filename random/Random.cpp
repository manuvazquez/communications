#include "Random.h"
#include <math.h>

using std::complex;

double Random::randn ()
{
  if(_havesmpl)
  {
	  _havesmpl = false;
	  return _smpl;
  }
  else
  {
	  double x1, x2, w;

	  do
		{
		  x1 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
		  x2 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
		  w = x1 * x1 + x2 * x2;
		}
	  while (w >= 1.0);

	  w = sqrt (-2.0 * log (w) / w);

	  _smpl = x2 * w;
	  _havesmpl = true;

	  return x1 * w;
  }
}

complex<double> Random::complexRandn()
{
  double x1, x2, w, y1, y2;

  do
    {
      x1 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
      x2 = 2.0 * (rand_r (&_seed) / (double) RAND_MAX) - 1.0;
      w = x1 * x1 + x2 * x2;
    }
  while (w >= 1.0);

  w = sqrt (-2.0 * log (w) / w);
  y1 = x1 * w;
  y2 = x2 * w;

  return complex<double> (y1, y2);
}

void Random::toOctaveFileStream(const std::vector <Random> &vector,std::string name,std::ofstream &f)
{
    f << "# name: " << name << std::endl <<"# type: struct" << std::endl << "# length: 3" << std::endl;
	
	toOctaveFileStreamAux3D(std::vector<std::vector<std::vector<Random> > >(1,std::vector<std::vector<Random> >(1,vector)),name,f);
}

void Random::toOctaveFileStream(const std::vector<std::vector <Random> >&matrix,std::string name,std::ofstream &f)
{
    f << "# name: " << name << std::endl <<"# type: struct" << std::endl << "# length: 3" << std::endl;
	
	toOctaveFileStreamAux3D(std::vector<std::vector<std::vector<Random> > >(1,matrix),name,f);
}

void Random::toOctaveFileStream(const std::vector<std::vector<std::vector <Random> > >&matrix,std::string name,std::ofstream &f)
{
	f << "# name: " << name << std::endl <<"# type: struct" << std::endl << "# length: 3" << std::endl;
	
	toOctaveFileStreamAux3D(matrix,name,f);
}

void Random::toOctaveFileStreamAux3D(const std::vector<std::vector<std::vector <Random> > >&matrix,std::string name,std::ofstream &f)
{
	std::stringstream dimensionsString;
	
	if(matrix.size()==1)
		dimensionsString << std::string("# rows: ") <<  matrix[0].size() << std::endl << std::string("# columns: ") << matrix[0][0].size() << std::endl;
	else
		dimensionsString << std::string("# ndims: 3") << std::endl << matrix[0].size() << " " << matrix[0][0].size() << " " << matrix.size() << std::endl;
	
	f << "# name: seed" << std::endl << "# type: cell" << std::endl << dimensionsString.str();
	for(uint k=0;k<matrix.size();k++)
		for(uint j=0;j<matrix[k][0].size();j++)
			for(uint i=0;i<matrix[k].size();i++)
				f << "# name: <cell-element>" << std::endl << "# type: scalar" << std::endl << matrix[k][i][j]._seed << std::endl;
	
	f << "# name: havesmpl" << std::endl << "# type: cell" << std::endl << dimensionsString.str();
	for(uint k=0;k<matrix.size();k++)
		for(uint j=0;j<matrix[k][0].size();j++)
			for(uint i=0;i<matrix[k].size();i++)
				f << "# name: <cell-element>" << std::endl << "# type: scalar" << std::endl << matrix[k][i][j]._havesmpl << std::endl;
	
	f << "# name: smpl" << std::endl <<"# type: cell" << std::endl << dimensionsString.str();
	for(uint k=0;k<matrix.size();k++)
		for(uint j=0;j<matrix[k][0].size();j++)
			for(uint i=0;i<matrix[k].size();i++)
				f << "# name: <cell-element>" << std::endl << "# type: scalar" << std::endl << matrix[k][i][j]._smpl << std::endl;
}
