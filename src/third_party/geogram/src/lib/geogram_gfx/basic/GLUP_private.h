/*
 *  Copyright (c) 2012-2016, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef __GEOGRAM_GFX_BASIC_GLUP_PRIVATE__
#define __GEOGRAM_GFX_BASIC_GLUP_PRIVATE__

#include <geogram_gfx/basic/GLUP.h>
#include <geogram_gfx/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h> // For cell descriptors / marching cells.
/**
 * \file geogram_gfx/basic/GLUP_private.h
 * \brief Implementation classes for the GLUP API.
 */

namespace GLUP {
    using namespace GEO;

    /**
     * \brief Computes the inverse of a 4x4 matrix.
     * \param[out] inv the computed inverse of \p m
     * \param[in] m pointer to the input matrix
     * \retval GL_TRUE if the matrix \p m is invertible
     * \retval GL_FALSE if the matrix \p m is singular. Then
     *   \p inv receives the transpose of the comatrix of \p m
     */
    GLboolean invert_matrix(GLfloat inv[16], const GLfloat m[16]);

    /**
     * \brief Computes the inverse of a 4x4 matrix.
     * \param[out] inv the computed inverse of \p m
     * \param[in] m pointer to the input matrix
     * \retval GL_TRUE if the matrix \p m is invertible
     * \retval GL_FALSE if the matrix \p m is singular. Then
     *   \p inv receives the transpose of the comatrix of \p m
     */
    GLboolean invert_matrix(GLdouble inv[16], const GLdouble m[16]);
    
    /**
     * \brief Computes the product of two 4x4 matrices
     * \param[out] out the computed product \p m1 * \p m2
     * \param[in] m1,m2 pointers to the input matrices
     */
    void mult_matrices(
        GLfloat out[16], const GLfloat m1[16], const GLfloat m2[16]
    );

    /**
     * \brief Computes the product of two 4x4 matrices
     * \param[out] out the computed product \p m1 * \p m2
     * \param[in] m1,m2 pointers to the input matrices
     */
    void mult_matrices(
        GLdouble out[16], const GLdouble m1[16], const GLdouble m2[16]
    );

    /**
     * \brief Computes the product of a 4x4 matrix and a vector.
     * \param[out] out the computed product \p m * \p v
     * \param[in] m pointer to the input matrix
     * \param[in] v pointer to the input vector
     */
    void mult_matrix_vector(
        GLfloat out[4], const GLfloat m[16], const GLfloat v[4]
    );


    /**
     * \brief Computes the product of the transpose of a 
     *   4x4 matrix and a vector.
     * \param[out] out the computed product \p m * \p v
     * \param[in] m pointer to the input matrix
     * \param[in] v pointer to the input vector
     * \TODO: it seems that in GLSL, w = M*v uses this one !!
     */
    void mult_transpose_matrix_vector(
        GLfloat out[4], const GLfloat m[16], const GLfloat v[4]
    );
    
    /**
     * \brief Computes the product of a 4x4 matrix and a vector.
     * \param[out] out the computed product \p m * \p v
     * \param[in] m pointer to the input matrix
     * \param[in] v pointer to the input vector
     */
    void mult_matrix_vector(
        GLdouble out[4], const GLdouble m[16], const GLdouble v[4] 
    );

    /**
     * \brief Transposes a matrix in-place.
     * \param[in,out] m a pointer to the 16 single-precision
     *  floating point coefficients of the matrix to be transposed.
     */
    void transpose_matrix(GLfloat m[16]);

    /**
     * \brief Transposes a matrix in-place.
     * \param[in,out] m a pointer to the 16 double-precision
     *  floating point coefficients of the matrix to be transposed.
     */
    void transpose_matrix(GLdouble m[16]);
    
    /**
     * \brief For debugging, outputs a matrix to the standard error.
     * \param[in] m the matrix to be displayed.
     */
    void show_matrix(const GLfloat m[16]);

    /**
     * \brief For debugging, outputs a vector to the standard error.
     * \param[in] v the vector to be displayed
     */
    void show_vector(const GLfloat v[4]);
    
    /**
     * \brief Resets a matrix to the identity matrix.
     * \param[out] out the matrix to be reset.
     */
    void load_identity_matrix(GLfloat out[16]);

    /**
     * \brief Copies a vector of floats.
     * \param[out] to a pointer to the destination vector
     * \param[in] from a const pointer to the source vector
     * \param[in] dim the number of components to copy
     */
    inline void copy_vector(GLfloat* to, const GLfloat* from, index_t dim) {
        Memory::copy(to, from, sizeof(GLfloat)*dim);
    }

    /**
     * \brief Copies a vector of doubles to a vector of floats.
     * \param[out] to a pointer to the destination vector
     * \param[in] from a const pointer to the source vector
     * \param[in] dim the number of components to copy
     */
    inline void copy_vector(GLfloat* to, const GLdouble* from, index_t dim) {
        for(index_t i=0; i<dim; ++i) {
            to[i] = GLfloat(from[i]);
        }
    }

    /**
     * \brief Copies a vector of floats to a vector of doubles.
     * \param[out] to a pointer to the destination vector
     * \param[in] from a const pointer to the source vector
     * \param[in] dim the number of components to copy
     */
    inline void copy_vector(GLdouble* to, const GLfloat* from, index_t dim) {
        for(index_t i=0; i<dim; ++i) {
            to[i] = GLdouble(from[i]);
        }
    }
    
    /**
     * \brief Normalizes a vector.
     * \param[in,out] v a pointer to the 3 coordinates of the 3d 
     *  vector to be normalized.
     */
    inline void normalize_vector(GLfloat v[3]) {
        GLfloat s = 1.0f / ::sqrtf(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;            
    }


    /**
     * \brief Gives the number of vertices for each GLUP primitive.
     * \details The array is indexed by the GLUP primitive type
     *  (one of GLUP_POINTS, GLUP_LINES, ...)
     */
    extern index_t nb_vertices_per_primitive[];
    
    /**********************************************************************/
    
    class Context;
    
    /**
     * \brief A Matrix stack.
     * \details There are threeo matrix stacks in a context,
     *  for modelview matrices, projection matrices and
     *  texture coordinates.
     */
    class MatrixStack {
    public:

        /**
         * \brief Maximum number of matrices in a stack.
         */
        static const int MAX_DEPTH=16;

        /**
         * \brief MatrixStack constructor.
         */
        MatrixStack() : top_(0) {
            load_identity_matrix(top());
        }

        /**
         * \brief Gets the matrix on the top of
         *  the stack.
         * \return a pointer to the coefficients of
         *  the matrix.
         */
        GLfloat* top() {
            return stack_[top_].data();
        }

        /**
         * \brief Pushes a copy of the top matrix.
         */
        void push() {
            geo_assert(top_ != MAX_DEPTH-1);
            GLfloat* from = top();
            ++top_;
            GLfloat* to = top();
            copy_vector(to, from, 16);
        }

        /**
         * \brief Removes a matrix from the top of
         *  the stack.
         */
        void pop() {
            geo_assert(top_ != 0);
            --top_;
        }

    protected:
        struct Matrix {
            GLfloat coeff[16];
            GLfloat* data() {
                return &coeff[0];
            }
        };
        
    private:
        Matrix stack_[MAX_DEPTH];
        index_t top_;
    };
        

    /**
     * \brief Number of vertices/colors/tex_coords in a GLUP buffer used
     *  by immediate mode.
     */
    static const index_t IMMEDIATE_BUFFER_SIZE = 1024*1024/16;

    /**
     * \brief Index of an ImmediateBuffer in the ImmediateState.
     */
    enum GLUPattribute {
        GLUP_VERTEX_ATTRIBUTE    = 0,
        GLUP_COLOR_ATTRIBUTE     = 1,
        GLUP_TEX_COORD_ATTRIBUTE = 2
    };

    /**
     * \brief A buffer used by GLUP in immediate mode.
     */
    class ImmediateBuffer {
        
    public:

        /**
         * ImmediateBuffer constructor.
         */
        ImmediateBuffer() :
            data_(nil),
            dimension_(0),
            is_enabled_(false),
            VBO_(0) {
            initialize(4); // TODO
        }

        /**
         * ImmediateBuffer destructor.
         */
        ~ImmediateBuffer() {
            delete[] data_;
            if(VBO_ != 0) {
                glDeleteBuffers(1, &VBO_);
                VBO_ = 0;
            }
        }

        void initialize(index_t dim) {
            data_ = new GLfloat[dim * IMMEDIATE_BUFFER_SIZE];
            dimension_ = dim;
            is_enabled_ = true; 
        }

        /**
         * \brief Enables this ImmediateBuffer.
         */
        void enable() {
            is_enabled_ = true;
        }

        /**
         * \brief Disables this ImmediateBuffer.
         */
        void disable() {
            is_enabled_ = false;
        }

        /**
         * \brief Tests whether this ImmediateBuffer is enabled.
         * \retval true if this ImmediateBuffer is enabled
         * \retval false otherwise
         */
        bool is_enabled() const {
            return is_enabled_;
        }
        
        /**
         * \brief Sets the current attribute value.
         * \param[in] x,y,z,w the component of the current attribute
         * \details Components past the dimension of the attribute
         *   are ignored (e.g., if dimension is 2, z and w are ignored).
         */
        void set_current(GLfloat x, GLfloat y, GLfloat z, GLfloat w) {
            current_[0] = x;
            current_[1] = y;
            current_[2] = z;
            current_[3] = w;
        }

        /**
         * \brief Copies the current attribute value to a specified
         *  vertex in this buffer.
         * \param[in] v the vertex index
         * \pre v < IMMEDIATE_BUFFER_SIZE
         */
        void copy_current_to(index_t v) {
            geo_debug_assert(v < IMMEDIATE_BUFFER_SIZE);
            if(is_enabled()) {
                copy_vector(element_ptr(v), current_, dimension());
            }
        }

        /**
         * \brief Copies this attribute from a vertex to another one.
         * \param[in] to index of the destination vertex
         * \param[in] from index of the source vertex
         * \pre to < IMMEDIATE_BUFFER_SIZE && from < IMMEDIATE_BUFFER_SIZE
         */
        void copy(index_t to, index_t from) {
            geo_debug_assert(to < IMMEDIATE_BUFFER_SIZE);
            geo_debug_assert(from < IMMEDIATE_BUFFER_SIZE);
            if(is_enabled()) {
                copy_vector(element_ptr(to), element_ptr(from), dimension());
            }
        }

        /**
         * \brief Gets the dimension of the attribute.
         * \return the number of components of the attribute
         */
        index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Gets the size of the memory used by the buffer.
         * \return the size of the buffer in bytes.
         */
        size_t size_in_bytes() const {
            return IMMEDIATE_BUFFER_SIZE * dimension() * sizeof(GLfloat);
        }
        
        /**
         * \brief Gets a pointer to one attribute value by index.
         * \param[in] v index of the vertex
         * \return a pointer to the attribute, i.e. an array of 
         *  \p dimension() GLfloats
         */
        GLfloat* element_ptr(index_t v) {
            geo_debug_assert(v < IMMEDIATE_BUFFER_SIZE);
            return data_ + v*dimension_;
        }

        /**
         * \brief Gets a pointer to the data.
         * \return a pointer to the first attribute. All the storage
         *  is contiguous in memory.
         */
        GLfloat* data() {
            return data_;
        }

        /**
         * \brief Gets the Vertex Buffer Object.
         * \return a modifiable reference to the Id of the Vertex Buffer
         *  Object. Can be zero if no VBO is used.
         */
        GLuint& VBO() {
            return VBO_;
        }
        
    private:
        GLfloat* data_;
        GLfloat current_[4];
        index_t dimension_; 
        bool is_enabled_;
        GLuint VBO_;
    };

    /**
     * \brief Stores all the buffers used to implement
     *  the immediate-mode interface.
     */
    class ImmediateState {
    public:
        /**
         * \brief ImmediateState constructor.
         */
        ImmediateState() :
            current_vertex_(0),
            max_current_vertex_(0),
            primitive_(GLUP_POINTS),
            VAO_(0) {
            
            buffer.resize(3);
            buffer[GLUP_VERTEX_ATTRIBUTE].initialize(4);
            buffer[GLUP_COLOR_ATTRIBUTE].initialize(4);
            buffer[GLUP_TEX_COORD_ATTRIBUTE].initialize(4);
            
            // Vertex is always enabled
            buffer[GLUP_VERTEX_ATTRIBUTE].enable();
        }

        /**
         * \brief ImmediateState destructor.
         */
        ~ImmediateState() {
            if(VAO_ != 0) {
                glDeleteVertexArrays(1, &VAO_);
                VAO_ = 0;
            }
        }
        
        /**
         * \brief Gets the Vertex Array Object.
         * \return a modifiable reference to the Id of the Vertex Array Object.
         *   Can be 0 if no VAO is used.
         */
        GLuint& VAO() {
            return VAO_;
        }
        
        /**
         * \brief Copies an element, i.e. all the attributes
         *  attached to a vertex.
         * \param[in] to index of the destination vertex
         * \param[in] from index of the source vertex
         * \details Only attributes that are enabled are copied.
         */
        void copy_element(index_t to, index_t from) {
            for(index_t i=0; i<buffer.size(); ++i) {
                buffer[i].copy(to, from);
            }
        }

        
        /**
         * \brief Configures the immediate state for rendering
         *  primitives of a given type.
         * \param[in] primitive type of the primitives to be rendered
         */
        void begin(GLUPprimitive primitive) {
            current_vertex_ = 0;
            max_current_vertex_ =
            IMMEDIATE_BUFFER_SIZE - (
                IMMEDIATE_BUFFER_SIZE %
                nb_vertices_per_primitive[primitive]
            );
            primitive_ = primitive;
        }

        /**
         * \brief Advances to the next vertex.
         * \details This copies all the current values of all enabled attributes
         *  to the current vertex position.
         */
        void next_vertex() {
            for(index_t i=0; i<buffer.size(); ++i) {
                buffer[i].copy_current_to(current_vertex_);
            }
            ++current_vertex_;
        }

        /**
         * \brief Tests whether the buffers are full.
         * \details When buffers are fulled, their contents need to be sent
         *  to OpenGL before calling reset(). These operations are done
         *  by the Context.
         */
        bool buffers_are_full() {
            return (current_vertex_ == max_current_vertex_);
        }

        /**
         * \brief Resets the current vertex index to zero.
         */
        void reset() {
            current_vertex_ = 0;
        }

        /**
         * \brief Gets the primitive currently rendered, i.e.
         *  the argument to the latest invocation of begin()
         */
        GLUPprimitive primitive() const {
            return primitive_;
        }

        /**
         * \brief Gets the number of vertices stored in the buffers.
         */
        index_t nb_vertices() const {
            return current_vertex_;
        }

        /**
         * \brief Gets the number of primitives stored in the buffers.
         */
        index_t nb_primitives() const {
            return current_vertex_ / nb_vertices_per_primitive[
                primitive_
            ];
        }
        
        vector<ImmediateBuffer> buffer;
        
    private:
        index_t current_vertex_;
        index_t max_current_vertex_;
        GLUPprimitive primitive_;
        GLuint VAO_;
    };
    

    /**********************************************************/
    
    /**
     * \brief Base class for representing GLUP state variables.
     */
    class StateVariableBase {
    public:

        /**
         * \brief StateVariableBase default constructor.
         */
        StateVariableBase() : address_(nil), context_(nil) {
        }

        /**
         * \brief StateVariableBase constructor.
         * \param[in] context a pointer to the GLUP Context
         * \param[in] name the name of the variable, without 
         *  "GLUPStateBlock." (it is prepended automatically).
         */
        StateVariableBase(
            Context* context, const char* name
        ) {
            initialize(context,name);
        }

        /**
         * \brief Initializes a StateVariableBase.
         * \param[in] context a pointer to the GLUP Context
         * \param[in] name the name of the variable, without 
         *  "GLUPStateBlock." (it is prepended automatically
         *  when searching for the variable in the state).
         */
        void initialize(Context* context, const char* name);

        /**
         * \brief Gets the name of this StateVariableBase.
         * \return a const reference to the name, without
         *  "GLUPStateBlock." prepended to it.
         */
        const std::string& name() const {
            return name_;
        }
        
    protected:
        friend class Context;
        
        /**
         * \brief Gets the address of the StateVariableBase.
         * \return a pointer to the variable in the client-side
         *  representation of the UBO.
         */
        Memory::pointer address() const {
            return address_;
        }
        
        /**
         * \brief Indicates that the variables in the context
         *  need to be sent to OpenGL.
         */
        void flag_uniform_buffer_as_dirty();

        Memory::pointer address_;
        Context* context_;
        std::string name_;
    };

    /**
     * \brief A GLUP state variable of a given type.
     * \tparam T the type of the state variable
     */
    template <class T> class StateVariable : public StateVariableBase {
    public:

        /**
         * \brief StateVariable default constructor.
         */
        StateVariable() {
        }

        /**
         * \brief StateVariableBase constructor.
         * \param[in] context a pointer to the GLUP Context
         * \param[in] name the name of the variable, without 
         *  "GLUPStateBlock." (it is prepended automatically).
         * \param[in] value initial value of the variable
         */
        StateVariable(
            Context* context, const char* name, T value
        ) : StateVariableBase(context, name) {
            set(value);
        }

        /**
         * \brief Initializes a StateVariable.
         * \param[in] context a pointer to the GLUP Context
         * \param[in] name the name of the variable, without 
         *  "GLUPStateBlock." (it is prepended automatically).
         * \param[in] value initial value of the variable
         */
        void initialize(Context* context, const char* name, T value) {
            StateVariableBase::initialize(context, name);
            set(value);
        }

        /**
         * \brief Gets the value.
         * \return the value of this StateVariable.
         */
        T get() const {
            return *reinterpret_cast<T*>(address_);
        }

        /**
         * \brief Sets the value.
         * \param[in] val the new value
         * \details flags the uniform buffer as dirty
         */
        void set(T val) {
            *reinterpret_cast<T*>(address_) = val;
            flag_uniform_buffer_as_dirty();
        }
    };

    /**
     * \brief A GLUP state variable that contains an array
     *  of floating points. This concerns both vectors and
     *  matrices.
     */
    class FloatsArrayStateVariable : public StateVariableBase {
    public:

        /**
         * \brief FloatsArrayStateVariable default constructor.
         */
        FloatsArrayStateVariable() {
        }

        /**
         * \brief FloatsArrayStateVariable constructor.
         * \param[in] context a pointer to the GLUP Context
         * \param[in] name the name of the variable, without 
         *  "GLUPStateBlock." (it is prepended automatically).
         */
        FloatsArrayStateVariable(
            Context* context, const char* name
        ) : StateVariableBase(context, name) {
        }

        /**
         * \brief Gets a pointer to the variable.
         * \return a const pointer to the first GLUPfloat stored
         *  in the variable
         */
        const GLUPfloat* get_pointer() const {
            return reinterpret_cast<const GLUPfloat*>(address_);
        }

        /**
         * \brief Gets a modifiable pointer to the variable.
         * \return a modifiable pointer to the first GLUPfloat stored
         *  in the variable
         * \details This flags the uniform buffer as dirty in the
         *  Context.
         */
        GLUPfloat* get_pointer() {
            // This is a non-const pointer, therefore it will be
            // probably modified by client code (else the 'const'
            // version of get_pointer() would have been called).
            flag_uniform_buffer_as_dirty();
            return reinterpret_cast<GLUPfloat*>(address_);
        }
    };

    /**
     * \brief A GLUP state variable that contains a vector. 
     * \details This corresponds to vec2, vec3, vec4 GLSL types.
     */
    class VectorStateVariable : public FloatsArrayStateVariable {
    public:

        /**
         * \brief VectorStateVariable default constructor.
         */
        VectorStateVariable() : dimension_(0) {
        }

        /**
         * \brief VectorStateVariable constructor.
         * \param[in] context a pointer to the GLUP Context
         * \param[in] name the name of the variable, without 
         *  "GLUPStateBlock." (it is prepended automatically)
         * \param[in] dimension 2 for vec2, 3 for vec3, 4 for vec4
         */
        VectorStateVariable(
            Context* context, const char* name, index_t dimension
        ) : FloatsArrayStateVariable(context, name), dimension_(dimension) {
            clear();
        }

        /**
         * \brief Initializes a VectorStateVariable.
         * \param[in] context a pointer to the GLUP Context
         * \param[in] name the name of the variable, without 
         *  "GLUPStateBlock." (it is prepended automatically)
         * \param[in] dimension 2 for vec2, 3 for vec3, 4 for vec4
         */
        void initialize(Context* context, const char* name, index_t dimension) {
            FloatsArrayStateVariable::initialize(context, name);
            dimension_ = dimension;
            clear();
        }

        /**
         * \brief Gets the dimension.
         * \return the number of components of this vector
         */
        index_t dimension() const {
            return dimension_;
        }
        
        /**
         * \brief Gets the value.
         * \param[out] x a pointer to an array of dimension()
         *  GLfloats, where to store the value
         */
        void get(GLUPfloat* x) const {
            Memory::copy(x, address_, sizeof(GLUPfloat)*dimension_);
        }

        /**
         * \brief Sets the value.
         * \param[in] x a const pointer to an array of dimension()
         *  GLfloats that contains the new value
         */
        void set(const GLUPfloat* x) {
            Memory::copy(address_, x, sizeof(GLUPfloat)*dimension_);
            flag_uniform_buffer_as_dirty();
        }

        /**
         * \brief clears the vector to its default value.
         * \details For vec2, default value is (0.0, 0.0), for
         *  vec3, it is (0.0, 0.0, 0.0) and for vec4 it is 
         *  (0.0, 0.0, 0.0, 1.0)
         */
        void clear() {
            Memory::clear(address_, sizeof(GLUPfloat)*dimension_);
            if(dimension_ == 4) {
                reinterpret_cast<GLUPfloat*>(address_)[3] = 1.0f;
            }
            flag_uniform_buffer_as_dirty();            
        }
            
    protected:
        index_t dimension_;
    };


    /**
     * \brief The set of state variables that represent GLUP uniform state.
     * \details This reflects the GLSL GLUP uniform state declaration.
     * \see Context::GLSL_uniform_state_declaration()
     */
    struct UniformState {
        vector< StateVariable<GLboolean> > toggle;
        vector< VectorStateVariable>       color;
        VectorStateVariable                light_vector;
        VectorStateVariable                light_half_vector;
        StateVariable<GLfloat>             mesh_width;
        StateVariable<GLfloat>             cells_shrink;
        StateVariable<GLint>               picking_mode;
        StateVariable<GLint>               picking_id;      
        StateVariable<GLint>               base_picking_id; 
        StateVariable<GLint>               clipping_mode;
        StateVariable<GLint>               texture_mode;
        StateVariable<GLint>               texture_type;        
        VectorStateVariable                clip_plane;
        VectorStateVariable                world_clip_plane;
        FloatsArrayStateVariable           modelview_matrix;
        FloatsArrayStateVariable           modelviewprojection_matrix;
        FloatsArrayStateVariable           normal_matrix;
        FloatsArrayStateVariable           texture_matrix;
    };
    
    /**********************************************************************/

    /**
     * \brief Stores the programs and vertex array object used to display
     *  a primitive of a given type.
     */
    struct PrimitiveInfo {

        static const index_t nb_toggles_configs = 33;
        
        /**
         * \brief PrimitiveInfo constructor.
         */
        PrimitiveInfo():
            GL_primitive(0),
            vertex_gather_mode_VAO(0),
            vertex_gather_mode(false),
            implemented(false) {
            for(index_t i=0; i<nb_toggles_configs; ++i) {
                program[i] = 0;
                program_initialized[i] = false;
            }
        }

        /**
         * \brief PrimitiveInfo destructor.
         * \details Deletes the programs and vertex array object if need be.
         */
        ~PrimitiveInfo() {
            for(index_t i=0; i<nb_toggles_configs; ++i) {            
                if(program[i] != 0) {
                    glDeleteProgram(program[i]);
                    program[i] = 0;
                }
            }
            if(vertex_gather_mode_VAO != 0) {
                glDeleteVertexArrays(1,&vertex_gather_mode_VAO);
                vertex_gather_mode_VAO = 0;
            }
        }
        
        GLenum GL_primitive;
        GLuint program[nb_toggles_configs];
        bool program_initialized[nb_toggles_configs];
        GLuint vertex_gather_mode_VAO;
        bool vertex_gather_mode;
        bool implemented;
    };
    
    /**********************************************************************/

    /**
     * \brief Implements the MarchingCells algorithm.
     * \details MarchingCell compute the intersection between
     *  a cell and a plane, using only combinatorial information.
     *  It uses the static tables from the Mesh class, so that cell
     *  numberings are coherent between storage and graphics.
     */
    class MarchingCell {
    public:

        /**
         * \brief MarchingCell constructor
         * \param[in] prim the GLUP volumetric primitive, should be one
         *  of GLUP_TETRAHEDRA, GLUP_HEXAHEDRA, GLUP_PRISMS, GLUP_PYRAMIDS.
         */
        MarchingCell(GLUPprimitive prim);

        /**
         * \brief MarchingCell destructor.
         */
        ~MarchingCell();

        /**
         * \brief Gets the number of vertices.
         * \return the number of vertices in a cell
         */
        index_t nb_vertices() const {
            return nb_vertices_;
        }

        /**
         * \brief Gets the number of edges.
         * \return the number of edges in a cell
         */
        index_t nb_edges() const {
            return nb_edges_;
        }


        /**
         * \brief Gets the number of configurations.
         * \return the number of configurations
         */
        index_t nb_configs() const {
            return nb_configs_;
        }
        
        /**
         * \brief Gets a vertex by edge index and local vertex index.
         * \param[in] e the index of the edge
         * \param[in] lv the local vertex index in the edge, one of 0,1
         * \return the vertex index 
         */
        index_t edge_vertex(index_t e, index_t lv) const {
            geo_debug_assert(e < nb_edges());
            geo_debug_assert(lv < 2);
            return edge_[e*2+lv];
        }

        /**
         * \brief Gets the number of intersected edges in a configuration.
         * \param[in] config the vertex configuration bitcode. The
         *  bit corresponding to vertex v is set if v is on the positive
         *  side of the intersection plane.
         * \return the number of intersected edges in configuration \p config
         */
        index_t config_size(index_t config) const {
            geo_debug_assert(config < nb_configs());
            return config_size_[config];
        }

        /**
         * \brief Gets the maximum configuration size.
         * \return the largest number of vertices in an intersection polygon
         */
        index_t max_config_size() const {
            return max_config_size_;
        }

        /**
         * \brief Gets the list of intersected edges in a configuration.
         * \param[in] config the vertex configuration bitcode. The
         *  bit corresponding to vertex v is set if v is on the positive
         *  side of the intersection plane.
         * \return a pointer to the array of edge indices that correspond to
         *  this configuration, of size config_size()
         */
        const index_t* config_edges(index_t config) const {
            geo_debug_assert(config < nb_configs());
            return config_ + config * nb_edges_;
        }

        /**
         * \brief Gets the GLSL declaration of marching cell uniform state.
         * \return a pointer to GLSL source code that declares this
         *  marching cell's uniform state.
         */
        const char* GLSL_uniform_state_declaration() const {
            return GLSL_uniform_state_declaration_.c_str();
        }

        /**
         * \brief Gets the GLSL declaration of the function that
         *  computes the intersections.
         * \return a pointer to GLSL source code.
         */
        const char* GLSL_compute_intersections() const {
            return GLSL_compute_intersections_.c_str();
        }

        
        /**
         * \brief Gets the binding point of the uniform buffer that
         *  contains the tables for the marching cell.
         * \return the uniform binding point
         */
        GLuint uniform_binding_point() const {
            return uniform_binding_point_;
        }

        /**
         * \brief Creates a Uniform Buffer Object that contains
         *  the tables for the marching cell.
         */
        GLuint create_UBO();

        /**
         * \brief Binds the uniform state marching cell variables
         *  to a given program.
         */
        void bind_uniform_state(GLuint program);
        
    protected:

        /**
         * \brief Computes the intersection polygon for a configuration.
         * \param[in] config the vertex configuration bitcode. The
         *  bit corresponding to vertex v is set if v is on the positive
         *  side of the intersection plane.
         */
        void compute_config(index_t config);
        
        /**
         * \brief Moves from a given halfedge to the next halfege.
         * \details The halfedge is refered to as a facet index and a
         *  local vertex index within the facet.
         * \param[in,out] f the index of the facet
         * \param[in,out] lv the local index of the vertex in the facet.
         */
        void move_to_next(index_t& f, index_t& lv) {
            lv = (lv+1) % desc_->nb_vertices_in_facet[f];
        }

        /**
         * \brief Gets the origin vertex of a halfedge.
         * \details The halfedge is refered to as a facet index and a
         *  local vertex index within the facet.
         * \param[in] f the index of the facet
         * \param[in] lv the local index of the vertex in the facet.
         * \return the index of the origin vertex of the halfedge
         */
        index_t origin_vertex(index_t f, index_t lv) {
            return desc_->facet_vertex[f][lv];
        }

        /**
         * \brief Gets the destination vertex of a halfedge.
         * \details The halfedge is refered to as a facet index and a
         *  local vertex index within the facet.
         * \param[in] f the index of the facet
         * \param[in] lv the local index of the vertex in the facet.
         * \return the index of the destination vertex of the halfedge
         */
        index_t destination_vertex(index_t f, index_t lv) {
            move_to_next(f, lv);
            return origin_vertex(f, lv);
        }

        /**
         * \brief Gets the edge index that corresponds to a given halfedge.
         * \details The halfedge is refered to as a facet index and a
         *  local vertex index within the facet.
         * \param[in] f the index of the facet
         * \param[in] lv the local index of the vertex in the facet.
         * \return the index of the edge.
         */
        index_t edge(index_t f, index_t lv) {
            index_t v1 = origin_vertex(f, lv);
            index_t v2 = destination_vertex(f, lv);
            index_t result = vv_to_e_[v1*desc_->nb_vertices+v2];
            geo_debug_assert(result != index_t(-1));
            return result;
        }

        /**
         * \brief Moves from a given halfedge to the opposite halfege.
         * \details The halfedge is refered to as a facet index and a
         *  local vertex index within the facet.
         * \param[in,out] f the index of the facet
         * \param[in,out] lv the local index of the vertex in the facet.
         */
        void move_to_opposite(index_t& f, index_t& lv);

        /**
         * \brief Tests whether a given edge is intersected.
         * \details The halfedge is refered to as a facet index and a
         *  local vertex index within the facet.
         * \param[in] f the index of the facet
         * \param[in] lv the local index of the vertex in the facet.
         * \param[in] config the vertex configuration bitcode. The
         *  bit corresponding to vertex v is set if v is on the positive
         *  side of the intersection plane.
         * \retval true if the edge is intersected
         * \retval false otherwise
         */
        bool edge_is_intersected(index_t f, index_t lv, index_t config) {
            index_t v1 = origin_vertex(f, lv);
            index_t v2 = destination_vertex(f, lv);
            bool v1_in = ((config & 1u<<v1) != 0);
            bool v2_in = ((config & 1u<<v2) != 0);
            return (v1_in != v2_in);
        }

        /**
         * \brief Tests whether a vertex configuration bitcode is
         *  ambiguous.
         * \details A configuration is ambiguous if it results in several
         *  intersection polygons.
         * \param[in] config the vertex configuration bitcode. The
         *  bit corresponding to vertex v is set if v is on the positive
         *  side of the intersection plane.
         * \retval true if the configuration is ambiguous
         * \retval false otherwise
         */
        bool config_is_ambiguous(index_t config);
        
        /**
         * \brief Gets the first intersected halfedge given a 
         *  vertex configuration.
         * \param[out] f the facet of the first intersected halfedge
         * \param[out] lv the local vertex index of the first intersected
         *  halfedge
         * \param[in] config the vertex configuration bitcode. The
         *  bit corresponding to vertex v is set if v is on the positive
         *  side of the intersection plane.
         * \retval true if there was an intersection
         * \retval false otherwise
         */
        bool get_first_edge(index_t& f, index_t& lv, index_t config);

    private:
        /**
         * \brief Forbids copy.
         */
        MarchingCell(const MarchingCell& rhs);
        
        /**
         * \brief Forbids copy.
         */
        MarchingCell& operator=(const MarchingCell& rhs);
        
    private:
        const CellDescriptor* desc_;
        index_t vv_to_e_[64];
        index_t nb_vertices_;
        index_t nb_configs_;
        index_t* config_size_;
        index_t max_config_size_;
        index_t* config_;
        index_t nb_edges_;
        index_t* edge_;
        std::string GLSL_uniform_state_declaration_;
        std::string GLSL_compute_intersections_;
        GLuint uniform_binding_point_;
        GLuint UBO_;
    };
    
    
    /**********************************************************************/
    
    /**
     * \brief GLUP context stores a Uniform Buffer Object with state
     *  variables similar to OpenGL's fixed functionality pipeline, and
     *  a set of Vertex Buffer Objects to emulate OpenGL's immediate mode.
     */
    class Context {
    public:
        /**
         * \brief Gets the GLSL declaration of GLUP uniform state.
         * \return a pointer to GLSL source code that declares 
         *  GLUP uniform state.
         * \details Can be used by client-code shaders that need to
         *  have access to the GLUP uniform state.
         */
        static const char* uniform_state_declaration();
        
        /**
         * \brief Context constructor.
         */
        Context();

        /**
         * \brief Context destructor.
         */
        virtual ~Context();


        /**
         * \brief Gets the profile name associated with this context.
         */
        virtual const char* profile_name() const = 0; 
        
        /**
         * \brief Tests whether a given GLUP primitive supports array mode.
         * \details If array mode is supported, then one can use glupDrawArray()
         *  and glupDrawElements() with the specified primitive.
         * \param[in] prim the primitive to be tested.
         * \retval true if array mode is supported with \p prim
         * \retval false otherwise
         */
        virtual bool primitive_supports_array_mode(GLUPprimitive prim) const;
        
        /**
         * \brief Creates the uniform state and GLSL programs.
         * \details This function may throw exceptions if GLSL 
         *  functionalities are not implemented in the OpenGL driver.
         */
        virtual void setup();

        
        /**
         * \brief Copies the relevant variables from OpenGL 
         *  fixed functionality pipeline to this GLUP context.
         * \param[in] which_attribute an '|'-combination of 
         *  GLUP_MATRICES_ATTRIBUTES_BIT, GLUP_COLORS_ATTRIBUTES_BIT,
         *  GLUP_LIGHTING_ATTRIBUTES_BIT, GLUP_CLIPPING_ATTRIBUTES_BIT or
         *  GLUP_ALL_ATTRIBUTES.
         */
        virtual void copy_from_GL_state(GLUPbitfield which_attributes);


        /**
         * \brief Copies the relevant variables from this GLUP
         *  context to OpenGL fixed functionality pipeline.
         * \param[in] which_attribute an '|'-combination of 
         *  GLUP_MATRICES_ATTRIBUTES_BIT, GLUP_COLORS_ATTRIBUTES_BIT,
         *  GLUP_LIGHTING_ATTRIBUTES_BIT, GLUP_CLIPPING_ATTRIBUTES_BIT or
         *  GLUP_ALL_ATTRIBUTES.
         */
        virtual void copy_to_GL_state(GLUPbitfield which_attributes);
        
        /**
         * \brief Binds GLUP uniform state to a program.
         * \param[in] program the id of the GLSL program
         * \details If the program uses GLUP, then it 
         *  binds the program to GLUP uniform state, else this
         *  function does nothing.
         */
        virtual void bind_uniform_state(GLuint program);

        /**
         * \brief Replaces the top of the current matrix stack
         *  with the specified matrix.
         * \param[in] m the matrix that will replace the top of
         *  the current matrix stack
         */
        void load_matrix(const GLfloat m[16]) {
            copy_vector(matrix_stack_[matrix_mode_].top(), m, 16);
            matrices_dirty_ = true;
        }

        /**
         * \brief Replaces the top of the current matrix stack
         *  with the identity matrix.
         */
        void load_identity() {
            load_identity_matrix(matrix_stack_[matrix_mode_].top());
            matrices_dirty_ = true;
        }

        /**
         * \brief Post-multiplies the top of the current matrix stack
         *   with the specified matrix.
         * \param[in] m the matrix that will post-multiply the
         *   top of the current matrix stack.
         * \see matrix_mode()
         */
        void mult_matrix(const GLfloat m[16]) {
            GLfloat product[16];
            mult_matrices(product,m,matrix_stack_[matrix_mode_].top());
            load_matrix(product);
        }

        /**
         * \brief Pushes a copy of the top of the current stack matrix
         *  onto the current stack matrix.
         * \see matrix_mode(), pop_matrix()
         */
        void push_matrix() {
            matrix_stack_[matrix_mode_].push();
        }

        /**
         * \brief Pops the top of the current stack matrix.
         */
        void pop_matrix() {
            matrix_stack_[matrix_mode_].pop();
            matrices_dirty_ = true;
        }
        
        /**
         * \brief Sets the current matrix stack.
         * \param[in] matrix one of GLUP_MODELVIEW, GLUP_PROJECT
         * \details This determines on which matrix stack set_matrix(),
         *  mult_matrix(), push_matrix() and pop_matrix() operate.
         */
        void set_matrix_mode(GLUPmatrix matrix) {
            matrix_mode_ = matrix;
        }

        /**
         * \brief Gets the current matrix stack.
         * \return The current matrix stack, i.e. 
         *  one of GLUP_MODELVIEW, GLUP_PROJECT
         */
        GLUPmatrix get_matrix_mode() const {
            return matrix_mode_;
        }

        /**
         * \brief Creates a new vertex in the immediate mode
         *  buffers.
         * \param[in] x,y,z,w the coordinates of the vertex
         * \details The color and texture coordinates of the new
         *  vertex are initialized from the current color and 
         *  current texture coordinates.
         */
        void immediate_vertex(
            GLfloat x, GLfloat y, GLfloat z=0.0f, GLfloat w=1.0f
        ) {
            immediate_state_.buffer[GLUP_VERTEX_ATTRIBUTE].set_current(x,y,z,w);
            immediate_state_.next_vertex();
            if(immediate_state_.buffers_are_full()) {
                flush_immediate_buffers();
            }
        }

        /**
         * \brief Specifies the current color for the immediate
         *  mode buffers.
         * \param[in] r,g,b,a the components of the current color.
         */
        void immediate_color(
            GLfloat r, GLfloat g, GLfloat b, GLfloat a = 1.0f
        ) {
            immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].set_current(r,g,b,a);
        }

        /**
         * \brief Specifies the current texture coordinates for the
         *  immediate mode buffers.
         * \param[in] s,t,u,v the current texture coordinates.
         */
        void immediate_tex_coord(
            GLfloat s, GLfloat t=0.0f, GLfloat u=0.0f, GLfloat v=1.0f
        ) {
            immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].set_current(
                s,t,u,v
            );
        }

        /**
         * \brief Sets the user program, to be used instead of
         *  the default GLUP programs for drawing the primitives.
         */
        void set_user_program(GLuint program) {
            user_program_ = program;
        }
        
        /**
         * \brief Begins rendering in immediate mode.
         * \param[in] primitive the primitive to be rendered.
         * \see immediate_vertex(), immediate_color(), immediate_tex_coord()
         */
        virtual void begin(GLUPprimitive primitive);

        /**
         * \brief Ends rendering in immediate mode.
         * \see begin()
         */
        virtual void end();

        /**
         * \brief Draws primitives using current OpenGL array bindings.
         * \details This function operates just like glDrawArrays(), 
         *  except that its \p primitive argument is a GLUPprimitive
         *  instead of regular OpenGL primitive. Internally it uses
         *  a (possibly different) OpenGL primitive, as well as
         *  a GLSL program to reinterpret it.
         * \param[in] primitive the GLUP primitive type
         * \param[in] first first index to be rendered
         * \param[in] count number of vertices to be rendered
         */
        virtual void draw_arrays(
            GLUPprimitive primitive, GLUPint first, GLUPsizei count
        );

        /**
         * \brief Draws primitives using current OpenGL array bindings.
         * \details This function operates just like glDrawElements(), 
         *  except that its \p primitive argument is a GLUPprimitive
         *  instead of regular OpenGL primitive. Internally it uses
         *  a (possibly different) OpenGL primitive, as well as
         *  a GLSL program to reinterpret it.
         * \param[in] primitive the GLUP primitive type
         * \param[in] count number of vertices to be rendered
         * \param[in] type type of element indices, as one of 
         *   GL_UNSIGNED_BYTE, GL_UNSIGNED_SHORT, or GL_UNSIGNED_INT
         * \param[in] indices a pointer to where the indices are stored.
         */
        virtual void draw_elements(
            GLUPprimitive primitive, GLUPsizei count,
            GLUPenum type, const GLUPvoid* indices
        );

        /**
         * \brief Gets a pointer to the representation of a uniform
         *  state variable in host memory from its (unqualified) name.
         * \param[in] name the name of the variable, without the suffix
         *   "GLUPStateBlock." (it is added internally by the function)
         * \return a pointer to where the variable is represented in
         *  client side.
         */
        virtual Memory::pointer get_state_variable_address(const char* name);

        /**
         * \brief Gets the uniform state.
         * \return a reference to the uniform state
         */
        UniformState& uniform_state() {
            return uniform_state_;
        }

        /**
         * \brief Gets the uniform state.
         * \return a const reference to the uniform state
         */
        const UniformState& uniform_state() const {
            return uniform_state_;
        }

        /**
         * \brief Indicates that the OpenGL representation
         *  of the uniform state is no longer in sync with
         *  the local copy.
         */
        void flag_uniform_buffer_as_dirty() {
            uniform_buffer_dirty_ = true;
        }


        /**
         * \brief Indicates that cached lighting information 
         *  needs to be recomputed.
         */
        void flag_lighting_as_dirty() {
            lighting_dirty_ = true;
        }


        /**
         * \brief Gets a pointer to the values of the matrix at the
         *  top of a given stack.
         * \param[in] matrix name of the stack, one of GLUP_MODELVIEW_MATRIX,
         *   GLUP_PROJECTION_MATRIX, GLUP_TEXTURE_MATRIX
         */
        GLUPfloat* get_matrix(GLUPmatrix matrix) {
            geo_debug_assert(matrix < 3);
            return matrix_stack_[matrix].top();
        }

    protected:

        /**
         * \brief Gets some GLSL declarations that depend on the current
         *  profile.
         * \return a const pointer to a string with the GLSL sources of
         *  the profile-dependent declarations.
         */
        const char* profile_dependent_declarations();

        /**
         * \brief This function is called before starting to
         *  render primitives. It is called by begin(), draw_arrays()
         *  and draw_elements().
         * \details Some primitives require to change some
         *  parameters in OpenGL. For instance, when we use
         *  GL_PATCH to gather the vertices of hexahedra and
         *  tetrahedra, the number of vertices per patch needs
         *  to be specified to OpenGL.
         */
        virtual void prepare_to_draw(GLUPprimitive primitive);
        
        /**
         * \brief Initializes the representation of the uniform state.
         */
        virtual void setup_state_variables();


        /**
         * \brief Set-ups the buffers for immediate rendering.
         * \details This creates VBOs and the VAO.
         */
        virtual void setup_immediate_buffers();

        /**
         * \brief Sends all the active immediate buffers to the GPU.
         * \details Overwrites the VBOs with the contents of the buffers.
         */
        virtual void stream_immediate_buffers();

        /**
         * \brief Setups the programs and VAOs used for each primitive.
         */
        virtual void setup_primitives();

        /**
         * \brief Setups GLSL programs for points.
         */
        virtual void setup_GLUP_POINTS();

        /**
         * \brief Setups GLSL programs for lines.
         */
        virtual void setup_GLUP_LINES();

        /**
         * \brief Setups GLSL programs for triangles.
         */
        virtual void setup_GLUP_TRIANGLES();

        /**
         * \brief Setups GLSL programs for quads.
         */
        virtual void setup_GLUP_QUADS();

        /**
         * \brief Setups GLSL programs for tetrahedra.
         */
        virtual void setup_GLUP_TETRAHEDRA();

        /**
         * \brief Setups GLSL programs for hexahedra.
         */
        virtual void setup_GLUP_HEXAHEDRA();

        /**
         * \brief Setups GLSL programs for prisms.
         */
        virtual void setup_GLUP_PRISMS();

        /**
         * \brief Setups GLSL programs for pyramids.
         */
        virtual void setup_GLUP_PYRAMIDS();

        /**
         * \brief Initializes the PrimitiveInfo associated with a 
         *  given GLUP primitive.
         * \param[in] glup_primitive the GLUP primitive.
         * \param[in] gl_primitive the GL primitive used by the implementation
         * \param[in] program the GLSL program used by the implementation
         */
        void set_primitive_info(
            GLUPprimitive glup_primitive, GLenum gl_primitive, GLuint program
        );

        /**
         * \brief Initializes the PrimitiveInfo associated with a 
         *  given GLUP primitive in vertex-gather mode.
         * \details In vertex-gather mode, all the coordinates of all vertices
         *  and all attributes of the primitive are gathered into a small
         *  number of vertices. This is required by 
         *  primitives that have a number of vertices that corresponds to no 
         *  existing OpenGL primitive (i.e., hexahedron and pyramid). 
         * \param[in] glup_primitive the GLUP primitive.
         * \param[in] gl_primitive the GL primitive used to display the GLUP
         *   primitive. The number of vertices of the GL primitive needs to
         *   be a divisor of the number of vertices of the GLUP primitive.
         * \param[in] program the GLSL program used by the implementation
         */
        void set_primitive_info_vertex_gather_mode(
            GLUPprimitive glup_primitive, GLenum gl_primitive, GLuint program
        );
        
        /**
         * \brief Flushes the immediate mode buffers.
         */
        virtual void flush_immediate_buffers();
            
        /**
         * \brief Copies GLUP uniform state to OpenGL.
         */
        void update_uniform_buffer();

        /**
         * \brief Updates the matrices in the uniform state
         *  from the matrices in the stacks.
         */
        void update_matrices();

        /**
         * \brief Updates the lighting in the uniform state.
         * \details Computes the half vector from the lighting
         *  vector.
         */
        void update_lighting();

        /**
         * \brief Updates the base picking id and sends it to
         *  OpenGL.
         */
        virtual void update_base_picking_id(GLint new_value);

        /**
         * \brief Gets the GLSL declaration of the functions that
         *   query the toggles from the state.
         * \details It can handle pre-determined fixed configurations
         *   of toggles depending on the arguments of 
         *   setup_shader_source_for_toggles(). In the initial configuration,
         *   all the toggles are read from the state.
         * \see setup_shaders_source_for_toggles()
         * \return a const char pointer to the GLSL declaration.
         */
        const char* toggles_declaration() const {
            return toggles_shader_source_.c_str();
        }

        /**
         * \brief Gets the GLSL declaration  of the constant that
         *  indicates the current primitive.
         * \return a string with the GLSL declaration.
         */
        std::string primitive_declaration(GLUPprimitive prim) const;
        
        /**
         * \brief Sets the string that describes the settings of
         *  the toggles for a given configuration.
         * \param[in] toggles_state an unsigned integer, with its bits
         *  corresponding to the state of each toggle
         * \param[in] toggles_undetermined an unsigned integer, with its bits
         *  set if the corresponding toggle state needs to be determined
         *  dynamically from GLUP state
         */
        void setup_shaders_source_for_toggles(
            GLUPbitfield toggles_state,
            GLUPbitfield toggles_undetermined=0
        );

        /**
         * \brief Sets the string that describes the settings of
         *  the toggles for a given configuration.
         * \param[in] toggles_config the identifier of the toggles
         *  configurations, used to index the GLSL program in the
         *  PrimitiveInfo class
         * \pre toggles_config < PrimitiveInfo::nb_toggles_configs
         */
        void setup_shaders_source_for_toggles_config(
            index_t toggles_config
        ) {
            if(toggles_config == (1 << GLUP_PICKING)) {
                geo_debug_assert(toggles_config == 32);
                setup_shaders_source_for_toggles(
                    (1 << GLUP_PICKING),  // picking=true
                    (1 << GLUP_CLIPPING)  // clipping=undecided (use state)
                );
            } else {
                setup_shaders_source_for_toggles(GLUPbitfield(toggles_config));
            }
        }
        
        /**
         * \brief Updates the toggles_config_ state variable from
         *  the individual state of each toggle.
         */
        void update_toggles_config();

        /**
         * \brief Creates the GLSL shader that corresponds to the
         *  specified primitive and current toggles configuration if
         *  not already initialized.
         * \param[in] primitive the primitive to be displayed
         */
        void create_program_if_needed(GLUPprimitive primitive);
        
    protected:
        
        // OpenGL Uniform state.
        GLuint default_program_;
        GLuint uniform_buffer_;
        GLuint uniform_binding_point_;
        GLint  uniform_buffer_size_;
        bool uniform_buffer_dirty_;

        // C++ Uniform state.
        Memory::byte* uniform_buffer_data_;        
        UniformState uniform_state_;

        bool lighting_dirty_;
        
        // Matrix stacks.
        GLUPmatrix matrix_mode_;
        MatrixStack matrix_stack_[3];
        bool matrices_dirty_;

        // Immediate mode buffers.
        ImmediateState immediate_state_;

        // Primitive informations (i.e., how to
        // draw a primitive of a given type).
        vector<PrimitiveInfo> primitive_info_;

        // The marching cells, for computing
        // intersections when clipping mode
        // is GLUP_CLIP_SLICE_CELLS
        MarchingCell marching_tet_;
        MarchingCell marching_hex_;
        MarchingCell marching_prism_;
        MarchingCell marching_pyramid_;
        
        GLuint user_program_;
        
        index_t toggles_config_;
        std::string toggles_shader_source_;

        bool precompile_shaders_;

        bool use_core_profile_;
        bool use_ES_profile_;
    };

    /*********************************************************************/

    /**
     * \brief Implementation of GLUP using modern OpenGL with GLSL 1.50
     *  shaders. 
     * \details All the primitives are implemented with good performance.
     *  Hexahedra and prisms do not support array mode (glupDrawArrays(),
     *  glupDrawElements()). This is because there is no standard OpenGL
     *  primitive with 8 or 5 vertices (except the configurable GL_PATCH
     *  that requires GLSL 4.40).
     */
    class Context_GLSL150 : public Context {
    public:
        /**
         * \copydoc Context::profile_name()
         */
        virtual const char* profile_name() const;

        /**
         * \copydoc Context::setup()
         */
        virtual void setup();
        
    protected:
        /**
         * \copydoc Context::setup_GLUP_POINTS()
         */
        virtual void setup_GLUP_POINTS();

        /**
         * \copydoc Context::setup_GLUP_LINES()
         */
        virtual void setup_GLUP_LINES();

        /**
         * \copydoc Context::setup_GLUP_TRIANGLES()
         */
        virtual void setup_GLUP_TRIANGLES();

        /**
         * \copydoc Context::setup_GLUP_QUADS()
         */
        virtual void setup_GLUP_QUADS();

        /**
         * \copydoc Context::setup_GLUP_TETRAHEDRA()
         */
        virtual void setup_GLUP_TETRAHEDRA();

        /**
         * \copydoc Context::setup_GLUP_PRISMS()
         */
        virtual void setup_GLUP_PRISMS();

        /**
         * \copydoc Context::setup_GLUP_HEXAHEDRA()
         */
        virtual void setup_GLUP_HEXAHEDRA();

        /**
         * \copydoc Context::setup_GLUP_PYRAMIDS()
         */
        virtual void setup_GLUP_PYRAMIDS();
    };

    /*********************************************************************/

    /**
     * \brief Implementation of GLUP using modern OpenGL with GLSL 4.40
     *  shaders. 
     * \details This mostly reuses the GLSL 1.50 implementation, except
     *  for hexahedra and prisms, where it uses a tessellation shader to
     *  fetch the vertices. This is because GL_PATCH has a configurable
     *  number of vertices.
     */
    class Context_GLSL440 : public Context_GLSL150 {
    public:
        /**
         * \copydoc Context::profile_name()
         */
        virtual const char* profile_name() const;
        
    protected:
        /**
         * \copydoc Context::setup_GLUP_HEXAHEDRA()
         */
        virtual void setup_GLUP_HEXAHEDRA();

        /**
         * \copydoc Context::setup_GLUP_PYRAMIDS()
         */
        virtual void setup_GLUP_PYRAMIDS();
    };

    /*********************************************************************/

    /**
     * \brief Implementation of GLUP using Vanilla (old-style) OpenGL.
     * \details This implementation does not use any shader. It is used
     *  as a fallback when the initialization of the other ones fails.
     *  Some primitive may be not implemented, degraded or of very low
     *  performance.
     */
    class Context_VanillaGL : public Context {
    public:
        /**
         * \brief Context_VanillaGL constructor.
         */
        Context_VanillaGL();

        /**
         * \copydoc Context::profile_name()
         */
        virtual const char* profile_name() const;
        
        /**
         * \copydoc Context::primitive_supports_array_mode()
         */
        virtual bool primitive_supports_array_mode(GLUPprimitive prim) const;

        /**
         * \copydoc Context::begin()
         */
        virtual void begin(GLUPprimitive primitive);

        /**
         * \copydoc Context::end()
         */
        virtual void end();
        
    protected:

        /**
         * \brief Configures texturing-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin()
         */
        void configure_OpenGL_texturing();

        /**
         * \brief Configures lighting-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin()
         */
        void configure_OpenGL_lighting();

        /**
         * \brief Configures clipping-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin(). 
         */
        void configure_OpenGL_clipping();
        
        /**
         * \brief Configures lighting-related OpenGL state
         *  variables according to the GLUP state variables.
         * \details This function is called by begin(). It needs
         *  to be called after configure_OpenGL_texturing() and
         *  configure_OpenGL_lighting() since it overrides texturing and
         *  lighting settings.
         */
        void configure_OpenGL_picking();        

        
        /**
         * \copydoc Context::setup()
         */
        virtual void setup();
        
        /**
         * \copydoc Context::flush_immediate_buffers()
         */
        virtual void flush_immediate_buffers();


        /**
         * \brief Flushes the immediate buffer with the
         *  current drawing modes. 
         * \details This function is separated from 
         *  flush_immediate_buffers(), since we need to
         *  flush the buffer twice when mesh drawing is
         *  enabled.
         */
        virtual void flush_immediate_buffers_once();
        
        /**
         * \copydoc Context::setup_immediate_buffers()
         */
        virtual void setup_immediate_buffers();


        /**
         * \copydoc Context::setup_primitives()
         */
        virtual void setup_primitives();
        
        /**
         * \copydoc Context::get_state_variable_address()
         */
        Memory::pointer get_state_variable_address(const char* name);


        /**
         * \brief Shrinks the cells in the immediate buffer.
         * \details Applies the shrinking factor (state variable
         *   "cells_shrink") to all the cells stored in the current
         *   immediate buffer. Since there is no function to query
         *   the content of the current buffer, modidying it is 
         *   acceptable.
         */
        void shrink_cells_in_immediate_buffers();

        /**
         * \brief Updates v_is_visible_[] according to
         *  current clipping plane.
         */
        void classify_vertices_in_immediate_buffers();

        /**
         * \brief Tests whether the cell starting at a given vertex
         *  in the immediate buffer is clipped, according to current
         *  clipping mode and current primitive type.
         * \param[in] first_v index of the first vertex of the cell in
         *  the immediate buffer
         * \retval true if the cell starting at \p first_v in the 
         *  immediate buffer is clipped-out
         * \retval false otherwise
         */
        bool cell_is_clipped(index_t first_v);

        /**
         * \brief Tests whether cells should be sliced.
         * \retval true if cells should be sliced
         * \retval false otherwise
         */
        bool clip_slice_cells() const {
            return (
                uniform_state_.toggle[GLUP_CLIPPING].get() &&
                uniform_state_.clipping_mode.get() == GLUP_CLIP_SLICE_CELLS &&
                immediate_state_.primitive() >= GLUP_TETRAHEDRA
            ) ;
        }
        
        /**
         * \brief Sends a vertex and its optional attributes to OpenGL.
         * \param[in] v the index of the vertex from the immediate buffer.
         */
        void output_vertex(index_t v) {
            if(immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].is_enabled()) {
                glColor4fv(
                    immediate_state_.buffer[GLUP_COLOR_ATTRIBUTE].element_ptr(v)
                );
            }
            if(immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].is_enabled()) {
                glTexCoord4fv(
                    immediate_state_.buffer[GLUP_TEX_COORD_ATTRIBUTE].
                    element_ptr(v)
                );
            }
            glVertex4fv(
                immediate_state_.buffer[GLUP_VERTEX_ATTRIBUTE].element_ptr(v)
            );
        }

        /**
         * \brief Sends a triangle normal to OpenGL
         * \param[in] v1,v2,v3 the indices of the three vertices from
         *  the immediate buffer.
         */
        void output_normal(index_t v1, index_t v2, index_t v3);

        /**
         * \brief Sends a quad normal to OpenGL
         * \param[in] v1,v2,v3,v4 the indices of the four vertices from
         *  the immediate buffer.
         */
        void output_normal(index_t v1, index_t v2, index_t v3, index_t v4);


        /**
         * \brief Sends a picking id to OpenGL and encodes it as a color.
         * \details The current base picking id is added to the id.
         *  If picking is deactivated or constant by object, 
         *  it does nothing.
         */
        void output_picking_id(index_t id) {
            if(pick_primitives_) {
                glPickingIdAsColor(
                    index_t(uniform_state_.base_picking_id.get()) + id
                );
            }
        }
        
        /**
         * \brief Sends a flat-shaded triangle to OpenGL
         * \param[in] v1,v2,v3 the indices of the three vertices from the
         *  immediate buffer.
         */
        void flat_shaded_triangle(index_t v1, index_t v2, index_t v3) {
            if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                output_normal(v1,v2,v3);
            }
            output_vertex(v1);
            output_vertex(v2);
            output_vertex(v3);
        }

        /**
         * \brief Sends a flat-shaded quad to OpenGL
         * \param[in] v1,v2,v3,v4 the indices of the three vertices from the
         *  immediate buffer.
         */
        void flat_shaded_quad(index_t v1, index_t v2, index_t v3, index_t v4) {
            if(uniform_state_.toggle[GLUP_LIGHTING].get()) {
                output_normal(v1,v2,v3,v4);
            }
            output_vertex(v1);
            output_vertex(v2);
            output_vertex(v4);
            output_vertex(v3);            
        }

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as point primitives.
         */
        void draw_immediate_buffer_GLUP_POINTS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as line primitives.
         */
        void draw_immediate_buffer_GLUP_LINES();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as triangle primitives.
         */
        void draw_immediate_buffer_GLUP_TRIANGLES();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as quad primitives.
         */
        void draw_immediate_buffer_GLUP_QUADS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as tetrahedra primitives.
         */
        void draw_immediate_buffer_GLUP_TETRAHEDRA();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as hexahedra primitives.
         */
        void draw_immediate_buffer_GLUP_HEXAHEDRA();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as prism primitives.
         */
        void draw_immediate_buffer_GLUP_PRISMS();

        /**
         * \brief Sends the contents of the immediate buffers to 
         *  OpenGL, as pyramid primitives.
         */
        void draw_immediate_buffer_GLUP_PYRAMIDS();

        /**
         * \brief Computes the intersection between the clipping plane and
         *  a segment.
         * \param[in] v1 index of the first extremity of the segment in the
         *  immediate buffer
         * \param[in] v2 index of the second extremity of the segment in the
         *  immediate buffer
         * \param[in] v1 index of where to wrote the intersection in the 
         *  isect_xxx arrays
         */
        void compute_intersection(index_t v1, index_t v2, index_t vi) {
            const GLUPfloat* eqn = world_clip_plane_;
            const GLUPfloat* p1 = immediate_state_.buffer[0].element_ptr(v1);
            const GLUPfloat* p2 = immediate_state_.buffer[0].element_ptr(v2);
            
            GLUPfloat t = -eqn[3] -(
                eqn[0]*p1[0] +
                eqn[1]*p1[1] +
                eqn[2]*p1[2]
            );

            GLUPfloat d =
                eqn[0]*(p2[0]-p1[0]) +
                eqn[1]*(p2[1]-p1[1]) +
                eqn[2]*(p2[2]-p1[2]) ;
            
            if(fabs(d) < 1e-6) {
                t = 0.5f;
            } else {
                t /= d;
            }

            GLUPfloat s = 1.0f - t;
            
            isect_point_[4*vi+0] = s*p1[0] + t*p2[0];
            isect_point_[4*vi+1] = s*p1[1] + t*p2[1];
            isect_point_[4*vi+2] = s*p1[2] + t*p2[2];
            isect_point_[4*vi+3] = 1.0f;
            
            if(immediate_state_.buffer[1].is_enabled()) {
                const GLUPfloat* c1 =
                    immediate_state_.buffer[1].element_ptr(v1);
                const GLUPfloat* c2 =
                    immediate_state_.buffer[1].element_ptr(v2);
                isect_color_[4*vi+0] = s*c1[0] + t*c2[0];
                isect_color_[4*vi+1] = s*c1[1] + t*c2[1];
                isect_color_[4*vi+2] = s*c1[2] + t*c2[2];
                isect_color_[4*vi+3] = s*c1[3] + t*c2[3];                
            }
            
            if(immediate_state_.buffer[2].is_enabled()) {
                const GLUPfloat* tex1 =
                    immediate_state_.buffer[2].element_ptr(v1);
                const GLUPfloat* tex2 =
                    immediate_state_.buffer[2].element_ptr(v2);
                
                isect_tex_coord_[4*vi+0] = s*tex1[0] + t*tex2[0];
                isect_tex_coord_[4*vi+1] = s*tex1[1] + t*tex2[1];
                isect_tex_coord_[4*vi+2] = s*tex1[2] + t*tex2[2];
                isect_tex_coord_[4*vi+3] = s*tex1[3] + t*tex2[3];
            }
        }

        /**
         * \brief Assemble the configuration code of a primitive
         *  relative to the clipping plane.
         * \param[in] first_v index of the first vertex of the 
         *  primitive in the immediate buffer
         * \param[in] nb_v number of vertices of the primitive
         * \return an integer with the i-th bit set if vertex i
         *  is visible, and unset if it is clipped.
         */
        index_t get_config(index_t first_v, index_t nb_v) {
            index_t result = 0;
            for(index_t lv=0; lv<nb_v; ++lv) {
                if(v_is_visible_[first_v+lv]) {
                    result = result | (1u << lv);
                }
            }
            return result;
        }

        /**
         * \brief Draws all the primitives from the immediate buffer using
         *  the marching cells algorithm.
         * \details This function is used when clipping is enabled and when
         *  clippping mode is GLUP_CLIP_SLICE_CELLS
         */
        void draw_immediate_buffer_with_marching_cells(
            const MarchingCell& cell
        );
        
    private:
        std::map<std::string, GLsizei> variable_to_offset_;

        /**
         * \brief Indicates for a given vertex whether it is clipped or
         *  is visible, according to the current clipping plane.
         */
        bool v_is_visible_[IMMEDIATE_BUFFER_SIZE];

        /**
         * \brief Indicates whether a picking id should be send to 
         *  OpenGL for each primitive.
         */
        bool pick_primitives_;

        /**
         * \brief computed intersections.
         * \details Used when clipping mode is GLUP_CLIP_SLICE_CELLS.
         */
        GLUPfloat isect_point_[12*4];

        /**
         * \brief computed colors of intersections.
         * \details Used when clipping mode is GLUP_CLIP_SLICE_CELLS.
         */
        GLUPfloat isect_color_[12*4];

        /**
         * \brief computed texture coordinates of intersections.
         * \details Used when clipping mode is GLUP_CLIP_SLICE_CELLS.
         */
        GLUPfloat isect_tex_coord_[12*4];


        /**
         * \brief Cached pointer to uniform state variable.
         */
        GLUPfloat* world_clip_plane_;
    };

    /*********************************************************************/

    
}

#endif
