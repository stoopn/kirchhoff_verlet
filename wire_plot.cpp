/* Copyright 2009 Norbert Stoop, stoopn@ethz.ch */

#include <string>
#include <iostream>
#include <fstream>
#include "wire.h"

void WireSim::output_vtk(std::string filename)
{
  std::ofstream posfile;
std::cout<<"h0\n";

  posfile.open(filename.c_str());
std::cout<<"h1\n";
  posfile<<"<?xml version=\"1.0\"?>\n";
  posfile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  posfile<<"<PolyData>\n";
  posfile<<"<Piece NumberOfPoints=\""<<total_number_of_nodes<<"\" NumberOfVerts=\"0\" NumberOfLines=\"1\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
  posfile<<"<Points>\n<DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (wireiter it=elems.begin(); it<elems.end(); it++)
    {
      posfile<<(float)(it->pos[0])<<" "<<(float)(it->pos[1])<<" "<<(float)(it->pos[2])<<"\n";
    }
std::cout<<"h2\n";
  posfile<<"</DataArray>\n</Points>\n";

  posfile<<"<Lines>\n<DataArray type=\"Int64\" Name=\"connectivity\"  format=\"ascii\">\n";
  for (int i=0; i<elems.size(); i++)
    {
      posfile<<i<<"\n";
    }
std::cout<<"h3\n";
  posfile<<"</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  posfile<<elems.size()<<"\n";
  posfile<<"</DataArray>\n</Lines>\n";
  posfile<<"<PointData>\n<DataArray type=\"Float64\" Name=\"E_bend\" format=\"ascii\">\n";
  for (int i=0; i<elems.size(); i++)
    {
      posfile<<(float)(elems[i].E_bend)<<"\n";
    }
std::cout<<"h4\n";
  posfile<<"</DataArray>\n";
  posfile<<"<DataArray type=\"Float64\" Name=\"E_tens\" format=\"ascii\">\n";
  for (int i=0; i<elems.size(); i++)
    {
      posfile<<(float)(elems[i].E_tens)<<"\n";
    }
std::cout<<"h5\n";
  posfile<<"</DataArray>\n";
  posfile<<"<DataArray type=\"Float64\" Name=\"E_tors\" format=\"ascii\">\n";
  for (int i=0; i<elems.size(); i++)
    {
      posfile<<(float)(elems[i].E_tors)<<"\n";
    }
  posfile<<"</DataArray>\n";
  posfile<<"<DataArray type=\"Float64\" Name=\"E_ct\" format=\"ascii\">\n";
  for (int i=0; i<elems.size(); i++)
    {
      posfile<<(float)(elems[i].E_ct)<<"\n";
    }
  posfile<<"</DataArray>\n</PointData>\n";
  posfile<<"</Piece>\n</PolyData>\n";
  posfile<<"</VTKFile>\n";
  posfile.close();
}
