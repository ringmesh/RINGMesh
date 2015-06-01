/*
 *  Copyright (c) 2012-2014, Bruno Levy
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

#include <geogram/basic/attributes.h>
#include <geogram/basic/permutation.h>
#include <algorithm>

namespace GEO {

    void AttributeStoreObserver::register_me(AttributeStore* store) {
        store->register_observer(this);
    }
    
    void AttributeStoreObserver::unregister_me(AttributeStore* store) {
        store->unregister_observer(this);        
    }

    
    AttributeStore::AttributeStore(
        index_t elemsize,
        index_t dim
    ) :
        element_size_(elemsize),
        dimension_(dim),
        cached_base_addr_(nil),
        cached_size_(0),
        lock_(0)
    {
    }
    
    void AttributeStore::notify(
        Memory::pointer base_addr, index_t size, index_t dim
    ) {
        if(
            size != cached_size_ ||
            base_addr != cached_base_addr_ ||
            dim != dimension_
        ) {
            cached_base_addr_ = base_addr;
            cached_size_ = size;
            dimension_ = dim;
            for(std::set<AttributeStoreObserver*>::iterator
                    it = observers_.begin(); it!=observers_.end(); ++it
            ) {
                (*it)->notify(cached_base_addr_, cached_size_, dim);
            }
        }
    }
    
    AttributeStore::~AttributeStore() {
        // It is illegal to keep an Attribute<> active
        // when the object it is bound to is destroyed.
        geo_assert(!has_observers());
    }

    void AttributeStore::register_observer(AttributeStoreObserver* observer) {
        Process::acquire_spinlock(lock_);
        geo_assert(observers_.find(observer) == observers_.end());
        observers_.insert(observer);
        observer->notify(cached_base_addr_, cached_size_, dimension_);
        Process::release_spinlock(lock_);        
    }

    void AttributeStore::unregister_observer(AttributeStoreObserver* observer) {
        Process::acquire_spinlock(lock_);        
        std::set<AttributeStoreObserver*>::iterator it =
            observers_.find(observer);
        geo_assert(it != observers_.end());
        observers_.erase(it);
        Process::release_spinlock(lock_);                
    }

    void AttributeStore::apply_permutation(
        const vector<index_t>& permutation
    ) {
        geo_debug_assert(permutation.size() <= cached_size_);
        Permutation::apply(
            cached_base_addr_, permutation, element_size_ * dimension_
        );
    }

    void AttributeStore::compress(
        const vector<index_t>& old2new
    ) {
        geo_debug_assert(old2new.size() <= cached_size_);
        index_t item_size = element_size_ * dimension_;
        for(index_t i=0; i<old2new.size(); ++i) {
            index_t j = old2new[i];
            if(j == index_t(-1) || j == i) {
                continue;
            }
            geo_debug_assert(j <= i);
            Memory::copy(
                cached_base_addr_+j*item_size,
                cached_base_addr_+i*item_size,
                item_size
            );
        }
    }

    void AttributeStore::zero() {
        Memory::clear(
            cached_base_addr_, element_size_ * dimension_ * cached_size_
        );
    }
    
    /*************************************************************************/

    AttributesManager::AttributesManager() : size_(0) {
    }

    AttributesManager::~AttributesManager() {
        clear(false,false);
    }
    
    void AttributesManager::resize(index_t new_size) {
        if(new_size == size_) {
            return;
        }
        for(
            std::map<std::string, AttributeStore*>::iterator
                it=attributes_.begin();
            it != attributes_.end(); ++it
        ) {
            it->second->resize(new_size);
        }
        size_ = new_size;
    }

    void AttributesManager::apply_permutation(
        const vector<index_t>& permutation
    ) {
        for(
            std::map<std::string, AttributeStore*>::iterator
                it=attributes_.begin();
            it != attributes_.end(); ++it
        ) {
            it->second->apply_permutation(permutation);
        }
    }

    void AttributesManager::compress(
        const vector<index_t>& old2new
    ) {
        for(
            std::map<std::string, AttributeStore*>::iterator
                it=attributes_.begin();
            it != attributes_.end(); ++it
        ) {
            it->second->compress(old2new);
        }
    }
    
    
    void AttributesManager::bind_attribute_store(
        const std::string& name, AttributeStore* as
    ) {
        geo_assert(find_attribute_store(name) == nil);
        attributes_[name] = as;
        as->resize(size_);
    }

    void AttributesManager::list_attribute_names(
        vector<std::string>& names
    ) const {
        names.clear();
        for(std::map<std::string, AttributeStore*>::const_iterator
                it = attributes_.begin(); it != attributes_.end();
            ++it) {
            names.push_back(it->first);
        }
    }
    
    AttributeStore* AttributesManager::find_attribute_store(
        const std::string& name
    ) {
        std::map<std::string, AttributeStore*>::iterator
            it = attributes_.find(name);
        if(it == attributes_.end()) {
            return nil;
        }
        return it->second;
    }


    void AttributesManager::delete_attribute_store(const std::string& name) {
        std::map<std::string, AttributeStore*>::iterator
            it = attributes_.find(name);
        geo_assert(it != attributes_.end());
        geo_assert(!it->second->has_observers());
        delete it->second;
        attributes_.erase(it);
    }

    void AttributesManager::delete_attribute_store(AttributeStore* as) {
        for(
            std::map<std::string, AttributeStore*>::iterator
                it=attributes_.begin();
            it != attributes_.end(); ++it
        ) {
            if(it->second == as) {
                delete as;
                attributes_.erase(it);
                return;
            }
        }
        geo_assert_not_reached;
    }


    void AttributesManager::clear(bool keep_attributes, bool keep_memory) {
        if(keep_attributes) {
            for(
                std::map<std::string, AttributeStore*>::iterator
                    it=attributes_.begin();
                it != attributes_.end(); ++it
            ) {
                it->second->clear(keep_memory);
            }
        } else {
            for(
                std::map<std::string, AttributeStore*>::iterator
                    it=attributes_.begin();
                it != attributes_.end(); ++it
            ) {
                delete it->second;
            }
            attributes_.clear();
        }
        size_ = 0;
    }

    void AttributesManager::zero() {
        for(
            std::map<std::string, AttributeStore*>::iterator
                it=attributes_.begin();
            it != attributes_.end(); ++it
        ) {
            it->second->zero();
        }
    }

    void AttributesManager::copy(const AttributesManager& rhs) {
        clear(false, false);
        resize(rhs.size());
        for(
            std::map<std::string, AttributeStore*>::const_iterator
                it=rhs.attributes_.begin();
            it != rhs.attributes_.end(); ++it
        ) {
            bind_attribute_store(it->first, it->second->clone());
        }
    }
    
    /************************************************************************/ 
    
}

