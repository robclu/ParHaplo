// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo link container 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_LINK_CONTAINER_CPU_HPP
#define PARHAPLO_LINK_CONTAINER_CPU_HPP

#include "link.hpp"

#include <utility>
#include <vector>

namespace haplo {
   
// ----------------------------------------------------------------------------------------------------------
/// @struct     LinkHasher
/// @brief      Hashing functor for links
/// @tparam     Type    The type used for the key pair
// ----------------------------------------------------------------------------------------------------------
template <typename Type>
struct LinkHasher {
    // -------------------------------------- ALIAS'S -------------------------------------------------------
    using key_type = std::pair<Type, Type>;
    // ------------------------------------------------------------------------------------------------------
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the hash of the key
    /// @param[in]  key     The key to hash
    // ------------------------------------------------------------------------------------------------------
    size_t operator()(const key_type key) const  
    {
        return std::hash<Type>()(key.first) ^ std::hash<Type>()(key.second);
    }
};

// Specialize for cpu's
template <>
class LinkContainer<devices::cpu> {
public:
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
    using key_type          = std::pair<size_t, size_t>;
    using link_map          = tbb::concurrent_unordered_map<key_type, Link, LinkHasher<size_t>>;
    using mapped_type       = typename link_map::mapped_type;
    using iterator          = typename link_map::iterator;
    using atomic_type       = tbb::atomic<size_t>;
    // ------------------------------------------------------------------------------------------------------
private:
    link_map      _links;                 //!< The links for the container
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor 
    // ------------------------------------------------------------------------------------------------------
    explicit LinkContainer() noexcept {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Insert a link into the link container 
    /// @param[in]  idx_one       The index of the first node
    /// @param[in]  idx_two       The index of the second node
    /// @param[in]  link            The link to insert into the container
    // ------------------------------------------------------------------------------------------------------
    inline void insert(const size_t idx_one, const size_t idx_two, const Link& link = Link()) 
    {
        _links[sort_pair(idx_one, idx_two)] = link;
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets if a key is present in the container 
    /// @param[in]  idx_two       The index of the lower node
    /// @param[in]  idx_two       The index of the upper node
    /// @return     If the link is in the map
    // ------------------------------------------------------------------------------------------------------    
    inline bool exists(const size_t idx_one, const size_t idx_two) const 
    {
        return (_links.find(sort_pair(idx_one, idx_two)) != _links.end());
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a link from the container -- no error checking
    /// @param[in]  idx_one       The index of the lower node
    /// @param[in]  idx_two       The index of the upper node
    /// @return     The link at the given positions, or an iterator to the end
    // ------------------------------------------------------------------------------------------------------    
    inline mapped_type& operator()(const size_t idx_one, const size_t idx_two) 
    {
        return _links[sort_pair(idx_one, idx_two)];
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a link from the container 
    /// @param[in]  idx_one       The index of the lower node
    /// @param[in]  idx_two       The index of the upper node
    /// @return     The link at the given positions, or an iterator to the end
    /// @throw      out_of_range    If the index is invalid
    // ------------------------------------------------------------------------------------------------------    
    inline mapped_type& at(const size_t idx_one, const size_t idx_two) 
    {
        return _links.at(sort_pair(idx_one, idx_two));
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a link from the container 
    /// @param[in]  idx_one       The index of the lower node
    /// @param[in]  idx_two       The index of the upper node
    /// @return     The link at the given positions, or an iterator to the end
    /// @throw      out_of_range    If the index is invalid
    // ------------------------------------------------------------------------------------------------------    
    inline const mapped_type& at(const size_t idx_one, const size_t idx_two) const
    {
        return _links.at(sort_pair(idx_one, idx_two));
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Takes two link indices for a links and sorts them into a pair
    /// @param[in]  idx_one     The first index
    /// @param[in]  idx_two     The second index
    // ------------------------------------------------------------------------------------------------------
    inline const std::pair<size_t, size_t> sort_pair(const size_t idx_one, const size_t idx_two) const 
    {
        return std::make_pair(std::min(idx_one, idx_two), std::max(idx_one, idx_two));
    }
};

}
#endif                  // PARAHAPLO_LINK_CONTAINER_CPU_HPP

