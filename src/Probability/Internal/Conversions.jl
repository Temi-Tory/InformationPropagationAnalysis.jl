# Helper function to convert from Float64 data to p-box data
function convert_to_pbox_data(
    node_priors::Dict{Int64, Float64},
    link_probability::Dict{Tuple{Int64, Int64}, Float64};
    uncertainty_type::Symbol = :none,  # Options: :none, :interval, :normal
    uncertainty_value::Float64 = 0.0
)
    # Convert node priors
    pbox_node_priors = Dict{Int64, pbox}()
    for (node, value) in node_priors
        if uncertainty_type == :interval && uncertainty_value > 0.0
            # Create interval p-box with fixed width uncertainty
            min_val = max(0.0, value - uncertainty_value)
            max_val = min(1.0, value + uncertainty_value)
            pbox_node_priors[node] = PBA.makepbox(PBA.interval(min_val, max_val))
        elseif uncertainty_type == :normal && uncertainty_value > 0.0
            # Create normal distribution with mean value and std of uncertainty_value
            pbox_node_priors[node] = PBA.normal(value, uncertainty_value)
            if PBA.minimum(pbox_node_priors[node]) < 0 || PBA.maximum(pbox_node_priors[node]) > 1
                left_bound = max(0.0, PBA.minimum(pbox_node_priors[node]))
                right_bound = min(1.0, PBA.maximum(pbox_node_priors[node]))
                pbox_node_priors[node] = PBA.makepbox(PBA.interval(left_bound, right_bound))
            end
        else
            # Create precise p-box (default)
            pbox_node_priors[node] = PBA.makepbox(PBA.interval(value, value))
        end
    end

    # Convert link probabilities
    pbox_link_probability = Dict{Tuple{Int64, Int64}, pbox}()
    for (edge, value) in link_probability
        if uncertainty_type == :interval && uncertainty_value > 0.0
            min_val = max(0.0, value - uncertainty_value)
            max_val = min(1.0, value + uncertainty_value)
            pbox_link_probability[edge] = PBA.makepbox(PBA.interval(min_val, max_val))
        elseif uncertainty_type == :normal && uncertainty_value > 0.0
            pbox_link_probability[edge] = PBA.normal(value, uncertainty_value)
            if PBA.minimum(pbox_link_probability[edge]) < 0 || PBA.maximum(pbox_link_probability[edge]) > 1
                left_bound = max(0.0, PBA.minimum(pbox_link_probability[edge]))
                right_bound = min(1.0, PBA.maximum(pbox_link_probability[edge]))
                pbox_link_probability[edge] = PBA.makepbox(PBA.interval(left_bound, right_bound))
            end
        else
            pbox_link_probability[edge] = PBA.makepbox(PBA.interval(value, value))
        end
    end

    return pbox_node_priors, pbox_link_probability
end
