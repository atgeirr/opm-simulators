/*
  Copyright 2026 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_SERIAL_COMMUNICATION_HEADER_INCLUDED
#define OPM_SERIAL_COMMUNICATION_HEADER_INCLUDED

#include <dune/common/parallel/communication.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm
{

class SerialCommunication
{
public:
    using Communicator = Dune::Communication<int>;

    class ParallelIndexSet
    {
    public:
        class LocalIndex
        {
        public:
            LocalIndex() = default;

            LocalIndex(std::size_t local,
                       int attribute = Dune::OwnerOverlapCopyAttributeSet::owner,
                       bool isPublic = true)
                : local_(local)
                , attribute_(attribute)
                , isPublic_(isPublic)
            {
            }

            std::size_t local() const
            {
                return local_;
            }

            int attribute() const
            {
                return attribute_;
            }

            bool isPublic() const
            {
                return isPublic_;
            }

        private:
            std::size_t local_ = 0;
            int attribute_ = Dune::OwnerOverlapCopyAttributeSet::owner;
            bool isPublic_ = true;
        };

        class Entry
        {
        public:
            Entry() = default;

            Entry(int global, LocalIndex local)
                : global_(global)
                , local_(std::move(local))
            {
            }

            int global() const
            {
                return global_;
            }

            const LocalIndex& local() const
            {
                return local_;
            }

        private:
            int global_ = 0;
            LocalIndex local_{};
        };

        using const_iterator = std::vector<Entry>::const_iterator;

        void beginResize()
        {
            entries_.clear();
            ++seqNo_;
        }

        void add(int global, const LocalIndex& local)
        {
            entries_.emplace_back(global, local);
        }

        void endResize()
        {
        }

        std::size_t size() const
        {
            return entries_.size();
        }

        int seqNo() const
        {
            return seqNo_;
        }

        const_iterator begin() const
        {
            return entries_.begin();
        }

        const_iterator end() const
        {
            return entries_.end();
        }

    private:
        std::vector<Entry> entries_;
        int seqNo_ = 0;
    };

    class RemoteIndices
    {
    public:
        using Neighbours = std::vector<int>;

        void setNeighbours(const Neighbours& neighbours)
        {
            neighbours_ = neighbours;
        }

        const Neighbours& getNeighbours() const
        {
            return neighbours_;
        }

        template <bool>
        void rebuild()
        {
        }

    private:
        Neighbours neighbours_;
    };

    class GlobalLookupIndexSet
    {
    public:
        GlobalLookupIndexSet() = default;

        explicit GlobalLookupIndexSet(const ParallelIndexSet& indexSet)
            : size_(indexSet.size())
        {
        }

        GlobalLookupIndexSet(const ParallelIndexSet&, std::size_t size)
            : size_(size)
        {
        }

        std::size_t size() const
        {
            return size_;
        }

    private:
        std::size_t size_ = 0;
    };

    explicit SerialCommunication(Dune::SolverCategory::Category category = Dune::SolverCategory::overlapping)
        : category_(category)
    {
    }

    explicit SerialCommunication(const Communicator& communicator,
                                 Dune::SolverCategory::Category category = Dune::SolverCategory::overlapping,
                                 bool /*freecomm*/ = false)
        : communicator_(communicator)
        , category_(category)
    {
    }

    Dune::SolverCategory::Category category() const
    {
        return category_;
    }

    const Communicator& communicator() const
    {
        return communicator_;
    }

    template <class T>
    void copyOwnerToAll(const T& source, T& dest) const
    {
        dest = source;
    }

    template <class T>
    void copyCopyToAll(const T& source, T& dest) const
    {
        dest = source;
    }

    template <class T>
    void addOwnerOverlapToAll(const T& source, T& dest) const
    {
        dest += source;
    }

    template <class T>
    void addOwnerCopyToAll(const T& source, T& dest) const
    {
        dest += source;
    }

    template <class T>
    void addOwnerCopyToOwnerCopy(const T& source, T& dest) const
    {
        dest += source;
    }

    template <class T1, class T2>
    void dot(const T1& x, const T1& y, T2& result) const
    {
        result = x.dot(y);
    }

    template <class T>
    auto norm(const T& x) const
    {
        return x.two_norm();
    }

    template <class T>
    void project(T&) const
    {
    }

    ParallelIndexSet& indexSet()
    {
        return indexSet_;
    }

    const ParallelIndexSet& indexSet() const
    {
        return indexSet_;
    }

    RemoteIndices& remoteIndices()
    {
        return remoteIndices_;
    }

    const RemoteIndices& remoteIndices() const
    {
        return remoteIndices_;
    }

    void buildGlobalLookup()
    {
        globalLookup_ = GlobalLookupIndexSet(indexSet_);
    }

    void buildGlobalLookup(std::size_t size)
    {
        globalLookup_ = GlobalLookupIndexSet(indexSet_, size);
    }

    void freeGlobalLookup()
    {
        globalLookup_ = GlobalLookupIndexSet{};
    }

    const GlobalLookupIndexSet& globalLookup() const
    {
        return globalLookup_;
    }

private:
    Communicator communicator_{};
    Dune::SolverCategory::Category category_;
    ParallelIndexSet indexSet_;
    RemoteIndices remoteIndices_;
    GlobalLookupIndexSet globalLookup_;
};

template <class Comm>
struct IsSerialCommunication : std::false_type
{
};

template <>
struct IsSerialCommunication<SerialCommunication> : std::true_type
{
};

template <class Comm>
inline constexpr bool is_serial_communication_v = IsSerialCommunication<std::remove_cv_t<std::remove_reference_t<Comm>>>::value;

} // namespace Opm

#endif // OPM_SERIAL_COMMUNICATION_HEADER_INCLUDED