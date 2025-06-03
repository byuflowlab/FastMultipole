# Efficient storage for N sets of M matrices (1x1, 2x2, ..., MxM per set)

struct MatrixSetStorage{T}
    data::Vector{T}      # All matrices, flattened and concatenated
    offsets::Vector{Int} # Starting index of each matrix in `data`
    sizes::Vector{Int}   # Size (n) for each matrix (n x n)
    N::Int               # Number of sets
    M::Int               # Max matrix size per set
end

function create_matrix_storage(N, M, T=Float64)
    total_matrices = N * M
    sizes = repeat(1:M, N)
    offsets = Vector{Int}(undef, total_matrices + 1)
    idx = 1
    for i in 1:total_matrices
        offsets[i] = idx
        n = sizes[i]
        idx += n * n
    end
    offsets[end] = idx
    data = Vector{T}(undef, offsets[end] - 1)
    MatrixSetStorage{T}(data, offsets, sizes, N, M)
end

# Access the k-th matrix of set i (1-based indices)
function get_matrix(storage::MatrixSetStorage, i::Int, k::Int)
    idx = (i - 1) * storage.M + k
    n = storage.sizes[idx]
    start = storage.offsets[idx]
    A = @view storage.data[start : start + n*n - 1]
    reshape(A, n, n)
end

# Example usage:
# storage = create_matrix_storage(N, M)
# mat = get_matrix(storage, set_index, matrix_index)

# ------------------------------------------------------------------------
# Efficiency comparison: KernelAbstractions vs. CUDA-customized kernels
#
# If you only use simple matrix operations (e.g., elementwise arithmetic, basic reductions),
# the performance difference between KernelAbstractions and a CUDA-customized kernel is usually small.
# KernelAbstractions can generate efficient code that is often close to CUDA.jl's performance,
# especially for memory-bound workloads.
#
# The main differences may appear for:
#   - Very small matrices (where kernel launch overhead dominates)
#   - Highly-tuned CUDA kernels using advanced features
#
# For most simple operations, expect performance within a few percent of CUDA.jl.
# Always benchmark for your specific use case.
