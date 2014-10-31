/*[
 * Association Scientifique pour la Geologie et ses Applications (ASGA)
 * Copyright (c) 1993-2013 ASGA. All Rights Reserved.
 *
 * This program is a Trade Secret of the ASGA and it is not to be:
 * - reproduced, published, or disclosed to other,
 * - distributed or displayed,
 * - used for purposes or on Sites other than described
 *   in the GOCAD Advancement Agreement,
 * without the prior written authorization of the ASGA. Licencee
 * agrees to attach or embed this Notice on all copies of the program,
 * including partial copies or modified versions thereof.
 ]*/

#ifndef __GRGMESH_SMART_POINTER__
#define __GRGMESH_SMART_POINTER__

#include <grgmeshlib/common.h>

namespace GRGMesh {

//____________________________________________________________________________

    /**
     * Automatic memory management using reference counting.
     * This class can be used with classes inheriting
     * the Counted class.
     * @see Counted
     */

    template< class T > class SmartPointer {

    public:
        SmartPointer()
            : pointer_( nil )
        {
        }
        SmartPointer( T* ptr )
            : pointer_( ptr )
        {
            T::ref( pointer_ ) ;
        }
        SmartPointer( const SmartPointer< T >& rhs )
            : pointer_( rhs )
        {
            T::ref( pointer_ ) ;
        }
        ~SmartPointer()
        {
            T::unref( pointer_ ) ;
        }

        SmartPointer< T >& operator=( T* ptr )
        {
            if( ptr != pointer_ ) {
                T::unref( pointer_ ) ;
                pointer_ = ptr ;
                T::ref( pointer_ ) ;
            }
            return *this ;
        }
        SmartPointer< T >& operator=( const SmartPointer< T >& rhs )
        {
            T* rhs_p = rhs ;
            if( rhs_p != pointer_ ) {
                T::unref( pointer_ ) ;
                pointer_ = rhs_p ;
                T::ref( pointer_ ) ;
            }
            return *this ;
        }

        /**
         * Makes the current instance a nil pointer ("forgets" the
         * current reference). 'p.forget();' is a shorthand for
         * 'p = nil;'
         */
        void forget()
        {
            T::unref( pointer_ ) ;
            pointer_ = nil ;
        }

        T* operator->() const
        {
            grgmesh_assert( pointer_ != nil ) ;
            return pointer_ ;
        }
        T& operator*() const
        {
            grgmesh_assert( pointer_ != nil ) ;
            return *pointer_ ;
        }
        operator T*() const
        {
            return pointer_ ;
        }

        /**
         * 'if(p.is_nil()) {...}' is a shorthand for 'if(p == nil) {...}'
         */
        bool is_nil() const
        {
            return ( pointer_ == nil ) ;
        }

    protected:
    private:
        T* pointer_ ;
    } ;

//____________________________________________________________________________

    /**
     * This is the base class to be used for objects having
     * "reference count" memory management. They can be
     * referred to by using SmartPointer<T>, calling ref()
     * and unref() when necessary.
     * @see SmartPointer
     */

    class GRGMESH_API Counted {

    public:
        Counted()
            : nb_refs_( 0 )
        {
        }
        virtual ~Counted()
        {
            grgmesh_assert( nb_refs_ == 0 ) ;
        }

        void ref() const
        {
            Counted* non_const_this = (Counted *) this ;
            non_const_this->nb_refs_++ ;
        }
        void unref() const
        {
            Counted* non_const_this = (Counted *) this ;
            non_const_this->nb_refs_-- ;

            grgmesh_assert( nb_refs_ >= 0 ) ;

            if( nb_refs_ == 0 ) {
                delete this ;
            }
        }
        bool is_shared() const
        {
            return ( nb_refs_ > 1 ) ;
        }

        static void ref( const Counted* counted )
        {
            if( counted != nil ) {
                counted->ref() ;
            }
        }
        static void unref( const Counted* counted )
        {
            if( counted != nil ) {
                counted->unref() ;
            }
        }

    private:
        int nb_refs_ ;
    } ;

}

#endif
