/*
 * Copyright (c) 2012-2018, Association Scientifique pour la Geologie et ses
 * Applications (ASGA). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of ASGA nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ASGA BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *     http://www.ring-team.org
 *
 *     RING Project
 *     Ecole Nationale Superieure de Geologie - GeoRessources
 *     2 Rue du Doyen Marcel Roubault - TSA 70605
 *     54518 VANDOEUVRE-LES-NANCY
 *     FRANCE
 */

#pragma once

#include <functional>
#include <future>

#include <ringmesh/basic/common.h>
#include <ringmesh/basic/types.h>

#include <geogram/basic/command_line.h>

/*!
 * @file  Future class for handling properly systems with no multithreading
 * @author Antoine Mazuyer
 */

namespace RINGMesh
{
    class basic_api TaskHandler
    {
        ringmesh_disable_copy_and_move( TaskHandler );

    public:
        TaskHandler() = default;

        TaskHandler( index_t nb_threads )
        {
            tasks_.reserve( nb_threads );
        }

        template < typename TASK, typename... Args >
        void execute( TASK&& task, const Args&... args )
        {
            auto to_execute = std::bind( std::forward< TASK&& >( task ),
                std::forward< const Args& >( args )... );
            if( multi_thread_ )
            {
                tasks_.emplace_back(
                    std::async( std::launch::async, to_execute ) );
            }
            else
            {
                to_execute();
            }
        }

        void wait_aysnc_tasks()
        {
            if( !tasks_.empty() )
            {
                for( auto& task : tasks_ )
                {
                    task.wait();
                }
            }
        }

    private:
        /// Vector containing all the tasks
        std::vector< std::future< void > > tasks_;

        /// Tells whether or not the multithreading
        /// is enabled.
        bool multi_thread_{ GEO::CmdLine::get_arg_bool( "sys:multithread" ) };
    };

    template < typename ACTION >
    void parallel_for( index_t size, const ACTION& action )
    {
        if( size == 0 )
        {
            return;
        }

        auto action_per_thread = [&action]( index_t start, index_t end ) {
            for( auto i : range( start, end ) )
            {
                action( i );
            }
        };

        index_t nb_threads{ std::min(
            size, std::thread::hardware_concurrency() ) };
        TaskHandler tasks{ nb_threads };
        index_t start{ 0 };

        index_t nb_tasks_per_thread{ size / nb_threads };
        for( auto thread : range( nb_threads - 1 ) )
        {
            ringmesh_unused( thread );
            tasks.execute(
                action_per_thread, start, start + nb_tasks_per_thread );
            start += nb_tasks_per_thread;
        }
        tasks.execute( action_per_thread, start, size );
        tasks.wait_aysnc_tasks();
    }

} // namespace RINGMesh
