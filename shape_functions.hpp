#ifndef ALIGN3D_SHAPE_FUNCTIONS__HPP
#define ALIGN3D_SHAPE_FUNCTIONS__HPP

#include "types.hpp"

namespace  Align3D {

  void Neighbor_Conformers(
			      const float* ref_coord,
			      const std::vector<double>& alpha_ref_vector,
			      const std::vector< unsigned >& ref_volumeAtomIndexVector,
			      const std::map< unsigned, std::vector< unsigned > >& ref_colorAtomType2IndexVectorMap,
			      const double ref_sov,
			      const double ref_sof,
			      const float* fit_coord,
			      const std::vector<double>& alpha_fit_vector,
			      const std::vector< unsigned >& fit_volumeAtomIndexVector,
			      const std::map< unsigned, std::vector< unsigned > >& fit_colorAtomType2IndexVectorMap,
			      const double fit_sov,
			      const double fit_sof,
			      const bool considerColorAtoms,
			      const bool is_qrat_nonzero,
			      const unsigned max_preiters,
			      const unsigned max_postiters,
			      const double opt_param,
			      float* matrix,
			      double& nbr_st,
			      double& nbr_ct
			      );


  // *** Gradients


  void ComputeOverlapAndAnalyticGradient(
					 double& o,   // overlap
					 double* g,   // gradient[ 7 ]
					 const double* m,   // rotation/translation matrix [ 12 ]
					 const float* ref_coord,    // per ref atom, XYZ coordinates [3*nref]
					 const std::vector<double>& alpha_ref_vector,
					 const std::vector< unsigned >& ref_volume_atom_index_vector,
					 const float* fit_coord,    // original per fit atom, XYZ coordinates [3*nfit]
					 const std::vector<double>& alpha_fit_vector,
					 const std::vector< unsigned >& fit_volume_atom_index_vector,
					 bool needGrad = true
					 );

  
  
  // ***  Interconvert quaternion/rotational matrix with translational matrix along for the ride
  
  void  VQuaternionTransArray2RotTransMatrix( double* m,          // double[ 12 ]
					      const double* q,    // double[ 7 ]
					      const bool normalize = false );
  void  VRotTransMatrix2QuaternionTransArray( double* q,          // double[ 7 ]
					      const double* m );  // double[ 12 ]

  
  void  VQuaternionTransArray2RotTransMatrix( float* m,          // float[ 12 ]
					      const double* q,    // float[ 7 ]
 					      const bool normalize = false );
  
  // *** Rotational/Translational matrix handling (3x3 matrix + 1x3 matrix == 12 floats)
       

  // Combine (matrix multiply) two rotational/translational matrices

 
  void  VCombineRotTransMatrixMatrix( double* t,           // double[ 12 ]
				      const double* t1,    // double[ 12 ]
				      const double* t2 );  // double[ 12 ]

  void  VRotTransMatrixNew( Uint8& packed,
			    const float* m );  // float[ 12 ]
  
  double ComputeFeatureOverlap( const float* ref_coord,                                                    // per ref atom, XYZ coordinates
				const std::vector<double>& alpha_ref_vector,
				const std::map< unsigned, std::vector< unsigned > >& ref_atom_type_to_index_vector_map,
				const float* fit_coord,                                                    // per fit atom, XYZ coordinates
				const std::vector<double>& alpha_fit_vector,
				const std::map< unsigned, std::vector< unsigned > >& fit_atom_type_to_index_vector_map,
				const double* m
				);

  double ComputeFeatureOverlap(
			       const float* ref_coord,
			       const std::vector<double>& alpha_ref_vector,
			       const std::map< unsigned, std::vector< unsigned > >& ref_atom_type_to_index_vector_map,
			       const float* fit_coord,
			       const std::vector<double>& alpha_fit_vector,
			       const std::map< unsigned, std::vector< unsigned > >& fit_atom_type_to_index_vector_map
			       );
  

  double ComputeShapeOverlap( const float* ref_coord,                        // per ref atom, XYZ coordinates
			      const std::vector<double>& alpha_ref_vector,
			      const std::vector< unsigned >& ref_volume_atom_index_vector,
			      const float* fit_coord,
			      const std::vector<double>& alpha_fit_vector,
			      const std::vector< unsigned >& fit_volume_atom_index_vector
			      );

  
  
  // Get Steric Center of Gravity
  template< class T >  void  Get3dStericCenterOfGravity( T* center,
							 const T* coord,
							 const unsigned ncoord )
  {
    if ( 0 == ncoord ) {
      center[ 0 ] = static_cast< T >( 0.0 );
      center[ 1 ] = static_cast< T >( 0.0 );
      center[ 2 ] = static_cast< T >( 0.0 );
      return;  // Nothing more to do...
    }

    // Loop over coordinates
    T  x = coord[ 0 ];
    T  y = coord[ 1 ];
    T  z = coord[ 2 ];
    for ( unsigned  icoord = 1;  icoord != ncoord;  ++icoord ) {
      x += coord[ icoord * 3 ];
      y += coord[ 1 + ( icoord * 3 ) ];
      z += coord[ 2 + ( icoord * 3 ) ];
    }

    center[ 0 ] = x / static_cast< T >( ncoord );
    center[ 1 ] = y / static_cast< T >( ncoord );
    center[ 2 ] = z / static_cast< T >( ncoord );
  }
  template< class T, class W >  void  Get3dCenterOfMass( T* center,
							 const T* coord,
							 const W* weight,
							 const unsigned ncoord )
  {
    if ( 0 == ncoord ) {
      center[ 0 ] = static_cast< T >( 0.0 );
      center[ 1 ] = static_cast< T >( 0.0 );
      center[ 2 ] = static_cast< T >( 0.0 );
      return;  // Nothing more to do...
    }

    // Loop over coordinates
    T  x = static_cast< T >( coord[ 0 ] * weight[ 0 ] );
    T  y = static_cast< T >( coord[ 1 ] * weight[ 0 ] );
    T  z = static_cast< T >( coord[ 2 ] * weight[ 0 ] );
    for ( unsigned  icoord = 1;  icoord != ncoord;  ++icoord ) {
      x += static_cast< T >( coord[     ( icoord * 3 ) ] * weight[ icoord ] );
      y += static_cast< T >( coord[ 1 + ( icoord * 3 ) ] * weight[ icoord ] );
      z += static_cast< T >( coord[ 2 + ( icoord * 3 ) ] * weight[ icoord ] );
    }

    center[ 0 ] = x / static_cast< T >( ncoord );
    center[ 1 ] = y / static_cast< T >( ncoord );
    center[ 2 ] = z / static_cast< T >( ncoord );
  }

  // Translate 3D Coordinates
  template< class T >  void  Translate3dCoords( T* coord,
						const T* translate,
						const unsigned ncoord )
  {
    T  x = translate[ 0 ];
    T  y = translate[ 1 ];
    T  z = translate[ 2 ];
    for( unsigned  icoord = 0;  icoord != ncoord;  ++icoord ) {
      coord[     ( icoord * 3 ) ] += x;
      coord[ 1 + ( icoord * 3 ) ] += y;
      coord[ 2 + ( icoord * 3 ) ] += z;
    }
  }

  // Translate 3D Coordinates to molecule center
  template< class T >  void  Translate3D2Center( T* coord,
						 const unsigned natoms )
  {
    // Get steric center
    T  cne[ 3 ];
    Get3dStericCenterOfGravity( cne, coord, natoms );

    // Translate to steric center
    cne[ 0 ] = -cne[ 0 ];
    cne[ 1 ] = -cne[ 1 ];
    cne[ 2 ] = -cne[ 2 ];
    Translate3dCoords( coord, cne, natoms );
  }



  
  void VStericFrameRotation( float* m,     // float[ 9 ]
			     const float* coord, // float[ n * 3 ]
			     const unsigned n,
			     float* ev );   // float[ 3 ]
 

  void VTransform2StericFrame( float* m,     // float[ 12 ]
			       float* coord, // float[ n * 3 ]
			       const unsigned n,
			       float* ev );   // float[ 3 ]


  void VApplyRotTransMatrix( float* c,
			     const unsigned n,
			     const float* m );

  // Avoids need to copy array
  void VApplyRotTransMatrix( float* c,
			     const float* d,
			     const unsigned n,
			     const float* m );

  

  void TransformCoordToStericFrameAndCalculateEigenvalues(
							     const unsigned natom,
							     float* coord,        // float[ natom * 3 ]
							     float* ev            // float[ 3 ]
							     );

  unsigned CalculateQrat( const float* ev, const float threshold );

  void TransformCoordToStericFrameAndCalculateQrat(
						      const unsigned natom,
						      const float threshold,
						      float* coord,        // float[ natom * 3 ]
						      unsigned& qrat
						      );
 
  void getJointColorTypeSet( const unsigned* ref_atom_types,
			     const unsigned ref_natom, 
			     const unsigned* fit_atom_types,
			     const unsigned fit_natom,
			     std::set<unsigned>& joint_color_atom_type_set );
    
    
  void getVolumeAtomIndexVector( const unsigned* atom_types,
				 const unsigned natom,
				 std::vector< unsigned >& atomIndexVector );
    
  void getColorAtomType2IndexVectorMap( const unsigned* atom_types,
					const unsigned natom,
					std::map< unsigned, std::vector< unsigned > >& atomType2IndexVectorMap );

  void restrictColorAtomType2IndexVectorMap( 
					    std::map< unsigned, std::vector< unsigned > >& atom_type_to_index_vector_map,
					    const std::set<unsigned>& joint_color_atom_type_set );
    
  void setUseCutOff( bool useCutOff );
  double getD2CutOff();
  double getAlpha( double r );
  double getA_ak( double alpha1, double alpha2 );
  double getOverlap0( double alpha1, double alpha2 );
  void setAlpha( const double* rad, const unsigned natom, std::vector< double >& alpha_vec );




    

}

#endif // ALIGN3D_SHAPE_FUNCTIONS__HPP

