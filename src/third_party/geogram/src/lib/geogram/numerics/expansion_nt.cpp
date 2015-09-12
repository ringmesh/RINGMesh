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

#include <geogram/numerics/expansion_nt.h>

namespace GEO {

    expansion_nt& expansion_nt::operator+= (const expansion_nt& rhs) {
        index_t e_capa = expansion::sum_capacity(rep(), rhs.rep());
        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_sum(rep(), rhs.rep());
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    expansion_nt& expansion_nt::operator+= (double rhs) {
        index_t e_capa = expansion::sum_capacity(rep(), rhs);

        // TODO: optimized in-place version to be used
        //   if(!shared() && e_capa < rep().capacity())

        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_sum(rep(), rhs);
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);

        return *this;
    }

    expansion_nt& expansion_nt::operator-= (const expansion_nt& rhs) {
        index_t e_capa = expansion::diff_capacity(rep(), rhs.rep());
        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_diff(rep(), rhs.rep());
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    expansion_nt& expansion_nt::operator-= (double rhs) {
        index_t e_capa = expansion::diff_capacity(rep(), rhs);

        // TODO: optimized in-place version to be used
        //   if(!shared() && e_capa < rep().capacity())

        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_diff(rep(), rhs);
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);

        return *this;
    }

    expansion_nt& expansion_nt::operator*= (const expansion_nt& rhs) {
        index_t e_capa = expansion::product_capacity(rep(), rhs.rep());
        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_product(rep(), rhs.rep());
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    expansion_nt& expansion_nt::operator*= (double rhs) {
        index_t e_capa = expansion::product_capacity(rep(), rhs);

        // TODO: optimized in-place version to be used
        //   if(!shared() && e_capa < rep().capacity())

        expansion* e = expansion::new_expansion_on_heap(e_capa);
        e->assign_product(rep(), rhs);
        expansion::unref_expansion(rep_);
        rep_ = e;
        expansion::ref_expansion(rep_);
        return *this;
    }

    /************************************************************************/

    expansion_nt expansion_nt::operator+ (const expansion_nt& rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::sum_capacity(rep(), rhs.rep())
        );
        e->assign_sum(rep(), rhs.rep());
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator- (const expansion_nt& rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::diff_capacity(rep(), rhs.rep())
        );
        e->assign_diff(rep(), rhs.rep());
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator* (const expansion_nt& rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::product_capacity(rep(), rhs.rep())
        );
        e->assign_product(rep(), rhs.rep());
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator+ (double rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::sum_capacity(rep(), rhs)
        );
        e->assign_sum(rep(), rhs);
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator- (double rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::diff_capacity(rep(), rhs)
        );
        e->assign_diff(rep(), rhs);
        return expansion_nt(e);
    }

    expansion_nt expansion_nt::operator* (double rhs) const {
        expansion* e = expansion::new_expansion_on_heap(
            expansion::product_capacity(rep(), rhs)
        );
        e->assign_product(rep(), rhs);
        return expansion_nt(e);
    }

    /************************************************************************/

    expansion_nt expansion_nt::operator- () const {
        expansion_nt result(*this);
        result.rep().negate();
        return result;
    }

    /************************************************************************/

    bool expansion_nt::operator> (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() > 0;
    }

    bool expansion_nt::operator>= (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() >= 0;
    }

    bool expansion_nt::operator< (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() < 0;
    }

    bool expansion_nt::operator<= (const expansion_nt& rhs) const {
        const expansion& diff = expansion_diff(rep(), rhs.rep());
        return diff.sign() <= 0;
    }

    bool expansion_nt::operator> (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() > 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() > 0;
    }

    bool expansion_nt::operator>= (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() >= 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() >= 0;
    }

    bool expansion_nt::operator< (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() < 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() < 0;
    }

    bool expansion_nt::operator<= (double rhs) const {
        if(rhs == 0.0) {
            return rep().sign() <= 0;
        }
        const expansion& diff = expansion_diff(rep(), rhs);
        return diff.sign() <= 0;
    }

    /************************************************************************/
}
