#include "xmlParser.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <iomanip> 
using namespace std;

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace boost;
using namespace boost::iostreams;

#ifndef __XML_BULDER_H__
#define __XML_BULDER_H__

struct Bond
    {

    Bond(unsigned int bond_type, unsigned int tag_a, unsigned int tag_b) : type(bond_type), a(tag_a), b(tag_b) { }
    unsigned int type;  
    unsigned int a;    
    unsigned int b;    
    };
	
struct Angle
    {

    Angle(unsigned int angle_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c) : type(angle_type), a(tag_a), b(tag_b), c(tag_c) { }
    unsigned int type;  //!< The type index of the angle
    unsigned int a;     //!< The tag of the first particle in the angle
    unsigned int b;     //!< The tag of the second particle in the angle
    unsigned int c;     //!< The tag of the third particle in the angle
    };

struct Dihedral
    {

    Dihedral(unsigned int dihedral_type, unsigned int tag_a, unsigned int tag_b, unsigned int tag_c, unsigned int tag_d) : type(dihedral_type), a(tag_a), b(tag_b), c(tag_c), d(tag_d) { }
    unsigned int type;  //!< The type index of the dihedral
    unsigned int a;     //!< The tag of the first particle in the dihedral
    unsigned int b;     //!< The tag of the second particle in the dihedral
    unsigned int c;     //!< The tag of the third particle in the dihedral
    unsigned int d;     //!< The tag of the forth particle in the dihedral
    };	

struct BoxDim
    {
    double xlo; //!< Minimum x coord of the box
    double xhi; //!< Maximum x coord of the box
    double ylo; //!< Minimum y coord of the box
    double yhi; //!< Maximum y coord of the box
    double zlo; //!< Minimum z coord of the box
    double zhi; //!< Maximum z coord of the box
    
    //! Constructs a useless box
    BoxDim();
    //! Constructs a box from -Len/2 to Len/2
    BoxDim(double Len);
    //! Constructs a box from -Len_x/2 to Len_x/2 for each dimension x
    BoxDim(double Len_x, double Len_y, double Len_z);
    };
struct vec
	{

		vec() : x(0.0), y(0.0), z(0.0)
		{
		}

		vec(double xp, double yp, double zp) : x(xp), y(yp), z(zp)
		{
		}
		double x;   //!< x-component
		double y;   //!< y-component
		double z;   //!< z-component		
	};
struct vec_int
	{

		vec_int() : x(0), y(0), z(0)
		{
		}

		vec_int(int xp, int yp, int zp) : x(xp), y(yp), z(zp)
		{
		}
		int x;   //!< x-component
		int y;   //!< y-component
		int z;   //!< z-component		
	};            

struct vec_uint
	{

		vec_uint() : x(0), y(0), z(0)
		{
		}

		vec_uint(unsigned int xp, unsigned int yp, unsigned int zp) : x(xp), y(yp), z(zp)
		{
		}
		unsigned int x;   //!< x-component
		unsigned int y;   //!< y-component
		unsigned int z;   //!< z-component		
	};	
const unsigned int NO_BODY = 0xffffffff;
const unsigned int NO_INDEX = 0xffffffff;


class Xmlbuilder 
    {
    public:
        //! Loads in the file and parses the data
        Xmlbuilder(const std::string &fname);

        //! Returns the number of particles to be initialized
        virtual unsigned int getNumDimensions() const;
        
        //! Returns the number of particles to be initialized
        virtual unsigned int getNumParticles() const;
        
        //! Returns the number of particle types to be initialized
        virtual unsigned int getNumParticleTypes() const;
        
        //! Returns the timestep of the simulation
        virtual unsigned int getTimeStep() const;
        
        //! Sets the timestep of the simulation
        virtual void setTimeStep(unsigned int ts);
        
        //! Returns the box the particles will sit in
        virtual BoxDim getBox() const;
        
        

        
        //! Returns the number of bond types to be created
        virtual unsigned int getNumBondTypes() const;
        
        //! Returns the number of angle types to be created
        virtual unsigned int getNumAngleTypes() const;
 
        //! Returns the number of dihedral types to be created
        virtual unsigned int getNumDihedralTypes() const;
        const std::vector< Bond >& getBond() const { return m_bonds; }
        const std::vector< Angle >& getAngle() const { return m_angles; }
        const std::vector< Dihedral >& getDihedral() const { return m_dihedrals; }
        //! Access the read particle positions
        const std::vector< vec >& getPos() const { return m_pos_array; }
	
         const std::vector< vec >& getVel() const { return m_vel_array; }       
        //! Access the read images
        const std::vector< vec_int >& getImage() const { return m_image_array; }
        const std::vector< double >& getMass() const { return m_mass_array; }			
        //! Access the read particle positions
        const std::vector< unsigned int >& getType() const { return m_type_array; }
		
        const std::vector< unsigned int >& getBody() const { return m_body_array; }		

        const std::vector< std::string >& getTypeMapping() const { return m_type_mapping; }	
        const std::vector< std::string >& getBondTypeMapping() const { return m_bond_type_mapping; }
        const std::vector< std::string >& getAngleTypeMapping() const { return m_angle_type_mapping; }
        const std::vector< std::string >& getDihedralTypeMapping() const { return m_dihedral_type_mapping; }		
        //! Access the read particle positions
        const std::vector< vec >& getOrientation() const { return m_orientation; }
        const std::string& getFilename() const { return m_fname; }		
    private:
        //! Helper function to read the input file
        void readFile(const std::string &fname);
        //! Helper function to parse the box node
        void parseBoxNode(const XMLNode& node);
        //! Helper function to parse the position node
        void parsePositionNode(const XMLNode& node);
        //! Helper function to parse the image node
        void parseImageNode(const XMLNode& node);
        //! Helper function to parse the velocity node
        void parseVelocityNode(const XMLNode& node);
        //! Helper function to parse the mass node
        void parseMassNode(const XMLNode& node);
        //! Helper function to parse diameter node
        void parseDiameterNode(const XMLNode& node);
        //! Helper function to parse the type node
        void parseTypeNode(const XMLNode& node);
        //! Helper function to parse the body node
        void parseBodyNode(const XMLNode& node);		
        //! Helper function to parse the bonds node
        void parseBondNode(const XMLNode& node);
        //! Helper function to parse the angle node
        void parseAngleNode(const XMLNode& node);
        //! Helper function to parse the dihedral node
        void parseDihedralNode(const XMLNode& node);
        //! Helper function to parse the improper node
        void parseMoleIdNode(const XMLNode& node);
        //! Parse charge node
        void parseChargeNode(const XMLNode& node);

        void parseOrientationNode(const XMLNode& node);
		
        //! Helper function for identifying the particle type id
        unsigned int getTypeId(const std::string& name);
        //! Helper function for identifying the bond type id
        unsigned int getBondTypeId(const std::string& name);
        //! Helper function for identifying the angle type id
        unsigned int getAngleTypeId(const std::string& name);
        //! Helper function for identifying the dihedral type id
        unsigned int getDihedralTypeId(const std::string& name);
        //! Helper function for identifying the improper type id
        unsigned int getMoleculerId(const std::string& name);
        
        std::map< std::string, boost::function< void (const XMLNode&) > > m_parser_map; //!< Map for dispatching parsers based on node type
        
        BoxDim m_box;   //!< Simulation box read from the file
        bool m_box_read;    //!< Stores the box we read in
        std::string m_fname;
        unsigned int m_num_dimensions;              //!< number of spatial dimensions
        std::vector< vec > m_pos_array;             //!< positions of all particles loaded
        std::vector< vec_int > m_image_array;       //!< images of all particles loaded
        std::vector< vec > m_vel_array;             //!< velocities of all particles loaded
        std::vector< double > m_mass_array;         //!< masses of all particles loaded
        std::vector< double > m_diameter_array;     //!< diameters of all particles loaded
        std::vector< unsigned int > m_type_array;   //!< type values for all particles loaded
        std::vector< unsigned int > m_body_array;   //!< body values for all particles loaded		
        std::vector< double > m_charge_array;       //!< charge of the particles loaded
        std::vector< Bond > m_bonds;                //!< Bonds read in from the file
        std::vector< Angle > m_angles;              //!< Angle read in from the file
        std::vector< Dihedral > m_dihedrals;        //!< Dihedral read in from the file
        unsigned int m_timestep;                    //!< The time stamp

        std::vector< vec > m_orientation;             //!< orientation of all particles loaded	
		
        std::vector<std::string> m_Molecular_mapping;
        std::vector<unsigned int> m_Molecular_id;         
        std::vector<std::string> m_type_mapping;          //!< The created mapping between particle types and ids
        std::vector<std::string> m_bond_type_mapping;     //!< The created mapping between bond types and ids
        std::vector<std::string> m_angle_type_mapping;    //!< The created mapping between angle types and ids
        std::vector<std::string> m_dihedral_type_mapping; //!< The created mapping between dihedral types and ids
        
    };



#endif



