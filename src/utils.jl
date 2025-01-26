struct DESY5SN_info
    data::DataFrame
    obs_flatdata::Vector{Float64}
    covariance::Matrix{Float64}
    precision::Matrix{Float64}
    std::Vector{Float64}

    function DESY5SN_info()

        DES_data = load_SN("DESY5")

        obs_flatdata = @. DES_data.data.MU - (5 * log10((1 + DES_data.data.zHEL) / (1 + DES_data.data.zHD)))
        covariance = DES_data.covariance + diagm(DES_data.data.MUERR_FINAL .^ 2)
        precision = inv(covariance)
        std = sqrt.(diag(covariance))

        new(DES_data.data, obs_flatdata, covariance, precision, std)
    end
end

struct PantheonPlusSN_info
    data::DataFrame
    obs_flatdata::Vector{Float64}
    covariance::Matrix{Float64}
    precision::Matrix{Float64}
    std::Vector{Float64}

    function PantheonPlusSN_info()

        PP_data = load_SN("PantheonPlus")

        z_mask = PP_data.data.zHD.>0.01
        masked_PP_df = PP_data.data[z_mask, :]

        obs_flatdata = @. masked_PP_df.m_b_corr - 5 * (log10((1 + masked_PP_df.zHEL) / (1 + masked_PP_df.zHD)))
        covariance = PP_data.covariance[z_mask, z_mask]
        precision = inv(covariance)
        std = sqrt.(diag(covariance))

        new(masked_PP_df, obs_flatdata, covariance, precision, std)
    end
end

struct Union3SN_info
    data::DataFrame
    obs_flatdata::Vector{Float64}
    covariance::Matrix{Float64}
    precision::Matrix{Float64}
    std::Vector{Float64}

    function Union3SN_info()

        un3_data = load_SN("Union3")

        obs_flatdata = un3_data.data.mb
        covariance = Hermitian(un3_data.covariance)
        precision = un3_data.inv_covariance
        std = sqrt.(diag(covariance))

        new(un3_data.data, obs_flatdata, covariance, precision, std)
    end
end
