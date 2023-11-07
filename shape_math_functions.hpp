#ifndef ALIGN3D_SHAPE_MATH_FUNCTIONS__HPP
#define ALIGN3D_SHAPE_MATH_FUNCTIONS__HPP


namespace  Align3D {

  template< class T > T expApproxDouble( T x ){
    const double C = 0.0009765625; // 1.0/1024; 1024 = 2^10;
    double x_double = static_cast< double >(x);
    double y = static_cast< double>(1) + C*x_double;
    y = y*y; // 1
    y = y*y; // 2
    y = y*y; // 3
    y = y*y; // 4
    y = y*y; // 5 
    y = y*y; // 6
    y = y*y; // 7
    y = y*y; // 8
    y = y*y; // 9
    y = y*y; // 10    
    return static_cast< T >( y );
  }
  
  template< class T > T expLookup( T x ){
    x = -500.0 * x;
    if( x <= (T) 0 ){
      return (T) 1;
    }else{
      long lookup = std::round( (T) x );
      if( lookup < expMinusSize ){
	return (T) expMinus[ lookup ];
      }
    }
    return (T) 0;
  }

  template< class T > T expLookupTrunc( T x ){
    x = -500.0 * x;
    if( x <= (T) 0 ){
      return (T) 1;
    }else{
      long lookup = std::round( (T) x );
      if( lookup < expMinusSizeTrunc ){
	return (T) expMinusTrunc[ lookup ];
      }
    }
    return (T) 0;
  }

  template< class T > T exp( T x ){
    return expApproxDouble( x );
    //    return expLookup( x );
    //    return expLookupTrunc( x );
    //    return std::exp( x );
    //    return (T) std::exp( (float) x );
  }
}

#endif // ALIGN3D_SHAPE_MATH_FUNCTIONS__HPP
  
