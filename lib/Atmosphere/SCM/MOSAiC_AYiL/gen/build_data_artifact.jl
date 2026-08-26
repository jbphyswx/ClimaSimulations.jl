"""
    build_data_artifact.jl

Fetch the MOSAiC AYiL model configuration and results from Zenodo, bind it as a
lazy artifact, and print the resulting `Artifacts.toml` entry.

Run once; the entry it prints is what ships in `Artifacts.toml`, after which the
package downloads the data on first use and nobody needs to run this again.

    > using Pkg: Pkg; Pkg.activate("gen"); include("gen/build_data_artifact.jl")

Zenodo record 10.5281/zenodo.10491362, CC-BY-4.0:
"A year in LES: Standardized daily high-resolution Large Eddy Simulations of the
Arctic boundary layer and clouds during the MOSAiC drift".
"""

using Pkg: Pkg
using SHA: SHA
using Downloads: Downloads

const NAME = "ayil_config_input_results"
const URL = "https://zenodo.org/api/records/10491362/files/ayil_config_input_results.zip/content"
# The record's published size. The sha256 computed below is what future downloads
# are verified against; this only catches a truncated fetch here.
const SIZE_BYTES = 910_740_303

artifacts_toml() = joinpath(dirname(@__DIR__), "Artifacts.toml")

"""
    main(; zip = joinpath(tempdir(), "\$NAME.zip"), keep_zip = true)

Download (unless `zip` already holds the right bytes), hash, unpack and bind.

The download is kept by default so that a failure later in the pipeline does not
cost another 911 MB.
"""
function main(; zip::AbstractString = joinpath(tempdir(), "$NAME.zip"), keep_zip::Bool = true)
    if isfile(zip) && filesize(zip) == SIZE_BYTES
        @info "reusing download" zip
    else
        @info "downloading" URL size_gb = round(SIZE_BYTES / 2^30, digits = 2) zip
        Downloads.download(URL, zip)
        got = filesize(zip)
        got == SIZE_BYTES ||
            error("downloaded $got bytes, the record says $SIZE_BYTES")
    end

    sha256 = bytes2hex(open(SHA.sha256, zip))
    @info "sha256" sha256

    # `create_artifact` hands us an empty directory and returns the tree hash of
    # whatever we leave in it. The record is a zip, so `Pkg`'s tar-based `unpack`
    # cannot read it.
    hash = Pkg.Artifacts.create_artifact() do dir
        @info "unpacking into" dir
        run(pipeline(`unzip -q -o $zip -d $dir`; stdout = devnull))
        # The zip carries a single top-level directory; lift the day directories
        # to the artifact root so the artifact *is* the data root.
        entries = readdir(dir)
        if length(entries) == 1 && isdir(joinpath(dir, entries[1]))
            inner = joinpath(dir, entries[1])
            for name in readdir(inner)
                mv(joinpath(inner, name), joinpath(dir, name))
            end
            rm(inner)
            @info "lifted" from = entries[1] days = length(readdir(dir))
        end
    end
    @info "artifact created" hash path = Pkg.Artifacts.artifact_path(hash)

    Pkg.Artifacts.bind_artifact!(
        artifacts_toml(),
        NAME,
        hash;
        download_info = [(URL, sha256)],
        lazy = true,
        force = true,
    )
    @info "bound into $(artifacts_toml())"
    keep_zip || rm(zip; force = true)
    println("\n----- Artifacts.toml -----")
    print(read(artifacts_toml(), String))
    return nothing
end

main()