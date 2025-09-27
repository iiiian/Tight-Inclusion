// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// High Resolution Timer.
//
// Resolution on Mac (clock tick)
// Resolution on Linux (1 us not tested)
// Resolution on Windows (clock tick not tested)

#pragma once

#include <tight_inclusion/config.hpp>

#include <chrono>
#include <cstddef>
#include <string>

////////////////////////////////////////////////////////////////////////////////

#ifdef TIGHT_INCLUSION_WITH_TIMER
#define TIGHT_INCLUSION_SCOPED_TIMER(total_time)                               \
    ticcd::ScopedTimer __tight_inclusion_timer(total_time)
#else
#define TIGHT_INCLUSION_SCOPED_TIMER(total_time)
#endif

namespace ticcd {

    class Timer {
    public:
        // default constructor
        Timer() : stopped(true), start_time(), end_time() {}
        // default destructor
        ~Timer() {}

        // start timer
        void start()
        {
            stopped = false; // reset stop flag
            start_time = clock::now();
        }

        // stop the timer
        void stop()
        {
            end_time = clock::now();
            stopped = true; // set timer stopped flag
        }
        // get elapsed time in second
        double getElapsedTime() { return this->getElapsedTimeInSec(); }
        // get elapsed time in second (same as getElapsedTime)
        double getElapsedTimeInSec()
        {
            return this->getElapsedTimeInMicroSec() * 0.000001;
        }

        // get elapsed time in milli-second
        double getElapsedTimeInMilliSec()
        {
            return this->getElapsedTimeInMicroSec() * 0.001;
        }
        // get elapsed time in micro-second
        double getElapsedTimeInMicroSec()
        {
            const auto now = stopped ? end_time : clock::now();
            const auto elapsed =
                std::chrono::duration<double, std::micro>(now - start_time);
            return elapsed.count();
        }

    private:
        // stop flag
        bool stopped;
        using clock = std::chrono::steady_clock;
        clock::time_point start_time;
        clock::time_point end_time;
    };

    class ScopedTimer {
    public:
        ScopedTimer() : m_total_time(nullptr) { start(); }

        ScopedTimer(double &total_time) : m_total_time(&total_time) { start(); }

        virtual ~ScopedTimer() { stop(); }

        inline void start() { m_timer.start(); }

        inline void stop()
        {
            m_timer.stop();
            if (m_total_time) {
                *m_total_time += getElapsedTimeInMicroSec();
            }
        }

        inline double getElapsedTimeInMicroSec()
        {
            return m_timer.getElapsedTimeInMicroSec();
        }

        inline const Timer &timer() { return m_timer; }

    protected:
        std::string m_msg;
        Timer m_timer;
        double *m_total_time;
    };

} // namespace ticcd
