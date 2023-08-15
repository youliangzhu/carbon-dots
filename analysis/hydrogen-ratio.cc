
#include "Xmlbuilder.h"

int main(int argc,char* argv[])
	{
	std::vector<std::string> filename;
	for(unsigned int i=0; i<(unsigned int)argc;i++ )
		{
		filename.push_back(argv[i]);
		}

	ofstream d("hydro-ratio.log");
	if (!d.good())
		{
		cerr << endl << "***Error! Unable to open dump file for writing: " << endl << endl;
		throw runtime_error("Error writting polymer_xml dump file");
		}		
	
	unsigned int nb_per_np = 538;
	unsigned int nframe=0;

	float deltr = 0.1;
	float rmax = 5.0;
	unsigned int ngmax = (unsigned int)(rmax/deltr);		
	std::vector<float> intra_dis;
	intra_dis.resize(ngmax);
	std::vector<float> inter_dis;
	inter_dis.resize(ngmax);
	std::vector<float> inter_dis3d;
	inter_dis3d.resize(ngmax*ngmax*ngmax);	
	std::vector<float> intra_dis3d;
	intra_dis3d.resize(ngmax*ngmax*ngmax);	
	unsigned int nball = 0; 
	for (unsigned int ii=1;ii<filename.size(); ii++)
		{
		string xml_open = filename[ii];
		Xmlbuilder build(xml_open.c_str());
		unsigned int N = build.getNumParticles();

		BoxDim box = build.getBox();
		std::vector<vec> pos = build.getPos();
		std::vector<vec> vel = build.getVel();
		std::vector<double> mass = build.getMass();
		std::vector<Bond> bonds = build.getBond();
		std::vector<Angle> angles = build.getAngle();		
		std::vector<unsigned int> type = build.getType();
		std::vector< std::string > typemapping = build.getTypeMapping();
		std::vector< std::string > bondtypemapping = build.getBondTypeMapping();
		std::vector< std::string > angletypemapping = build.getAngleTypeMapping();		
		std::vector<vec_int> image = build.getImage();

		float Lx = box.xhi - box.xlo;
		float Ly = box.yhi - box.ylo;
		float Lz = box.zhi - box.zlo;
		float Lxinv = 1.0/Lx;     
		float Lyinv = 1.0/Ly;
		float Lzinv = 1.0/Lz;

		unsigned int inter = 0;
		unsigned int intra = 0;
		
		for (unsigned int i = 0; i< bonds.size(); i++)
			{
			Bond bi = bonds[i];
			if(bondtypemapping[bi.type]=="sticky")
				{
				unsigned int a = (unsigned int)(bi.a/nb_per_np);
				unsigned int b = (unsigned int)(bi.b/nb_per_np);
				if (a == b)
					intra += 1;
				else
					inter += 1;				
				}
			}
			
		nball = N/nb_per_np;
		d<<nframe<<", inter, "<<float(inter)/float(nball)<<", intra, "<<float(intra)/float(nball)<<", inter-intra, "<<float(inter)/float(intra)<<endl;

		std::vector<vec> cent;
		cent.resize(nball);

		for (unsigned int i = 0; i< nball; i++)
			{
			vec cm = vec(0.0, 0.0, 0.0);
			vec p0 = pos[i*nb_per_np];
			
			for (unsigned int j = 0; j< 250; j++)
				{
				unsigned int idx = i*nb_per_np + j;
				float dx = pos[idx].x - p0.x;
				float dy = pos[idx].y - p0.y;
				float dz = pos[idx].z - p0.z;
				dx -= Lx*rint(dx*Lxinv);
				dy -= Ly*rint(dy*Lyinv);
				dz -= Lz*rint(dz*Lzinv);				
				cm.x += dx;
				cm.y += dy;
				cm.z += dz;				
				}
			cm.x /= float(250);
			cm.y /= float(250);
			cm.z /= float(250);	
			cm.x += p0.x;
			cm.y += p0.y;
			cm.z += p0.z;			
			cent[i] = cm; 
			}
			
			
		for (unsigned int i = 0; i< bonds.size(); i++)
			{
			Bond bi = bonds[i];
			if(bondtypemapping[bi.type]=="sticky")
				{
				unsigned int a = (unsigned int)(bi.a/nb_per_np);
				unsigned int b = (unsigned int)(bi.b/nb_per_np);
				vec pi = pos[bi.a];
				vec pj = pos[bi.b];
				vec cma = cent[a];
				vec cmb = cent[b];				
				float dx1 = pi.x - cma.x;
				float dy1 = pi.y - cma.y;
				float dz1 = pi.z - cma.z;
				dx1 -= Lx*rint(dx1*Lxinv);
				dy1 -= Ly*rint(dy1*Lyinv);
				dz1 -= Lz*rint(dz1*Lzinv);				
				float r1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
				float dx2 = pj.x - cmb.x;
				float dy2 = pj.y - cmb.y;
				float dz2 = pj.z - cmb.z;
				dx2 -= Lx*rint(dx2*Lxinv);
				dy2 -= Ly*rint(dy2*Lyinv);
				dz2 -= Lz*rint(dz2*Lzinv);				
				float r2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);				
				if (a == b)
					{
					if(r1<rmax&&abs(dx1)<rmax&&abs(dy1)<rmax&&abs(dz1)<rmax)
						{
						intra_dis[int(r1/deltr)] += 1;
						int ix1 = int(abs(dx1)/deltr);
						int iy1 = int(abs(dy1)/deltr);
						int iz1 = int(abs(dz1)/deltr);						
						intra_dis3d[ix1 + iy1*ngmax + iz1*ngmax*ngmax] += 1;
						}
					else
						cout<<r1<<" exceed rmax "<<endl;
					
					if(r2<rmax&&abs(dx2)<rmax&&abs(dy2)<rmax&&abs(dz2)<rmax)
						{
						intra_dis[int(r2/deltr)] += 1;
						int ix2 = int(abs(dx2)/deltr);
						int iy2 = int(abs(dy2)/deltr);
						int iz2 = int(abs(dz2)/deltr);						
						intra_dis3d[ix2 + iy2*ngmax + iz2*ngmax*ngmax] += 1;
						}					
					else
						cout<<r2<<" exceed rmax "<<endl;
					}
				else
					{
					if(r1<rmax&&abs(dx1)<rmax&&abs(dy1)<rmax&&abs(dz1)<rmax)
						{
						inter_dis[int(r1/deltr)] += 1;
						int ix1 = int(abs(dx1)/deltr);
						int iy1 = int(abs(dy1)/deltr);
						int iz1 = int(abs(dz1)/deltr);						
						inter_dis3d[ix1 + iy1*ngmax + iz1*ngmax*ngmax] += 1;
						}
					else
						cout<<r1<<" exceed rmax "<<endl;
					
					if(r2<rmax&&abs(dx2)<rmax&&abs(dy2)<rmax&&abs(dz2)<rmax)
						{
						inter_dis[int(r2/deltr)] += 1;
						int ix2 = int(abs(dx2)/deltr);
						int iy2 = int(abs(dy2)/deltr);
						int iz2 = int(abs(dz2)/deltr);						
						inter_dis3d[ix2 + iy2*ngmax + iz2*ngmax*ngmax] += 1;
						}					
					else
						cout<<r2<<" exceed rmax "<<endl;
					}			
				}
			}			
		nframe +=1;
		}
	
	ofstream to("hb_dis.log");
	for (unsigned int i =0; i<ngmax ; i++)
		{
		to <<float(i)*deltr <<" "<<intra_dis[i]/float(nframe*nball)<<"  "<<inter_dis[i]/float(nframe*nball)<<endl;
		}	

	ofstream to3("hb_dis3d.log");
	for (unsigned int i =0; i<ngmax ; i++)
		{
		for (unsigned int j =0; j<ngmax ; j++)
			{
			for (unsigned int k =0; k<ngmax ; k++)
				{			
				to3 <<float(i)*deltr <<" "<<float(j)*deltr <<" "<<float(k)*deltr <<" "<<intra_dis3d[i + j*ngmax + k*ngmax*ngmax]/float(nframe*nball)<<"  "<<inter_dis3d[i + j*ngmax + k*ngmax*ngmax]/float(nframe*nball)<<endl;
				}
			}
		}
	
	
	if (!d.good())
		{
		cerr << endl << "***Error! Unexpected error writing polymer dump file" << endl << endl;
		throw runtime_error("Error writting polymer dump file");
		}
	d.close();	
	 }

