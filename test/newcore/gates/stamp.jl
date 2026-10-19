# The commit a gate table is measured on, read when the script starts: the
# short hash of HEAD, and whether the working tree differed from it outside
# the tables the gate scripts write.
const COMMIT_STAMP = let repo = joinpath(@__DIR__, "..", "..", "..")
    head = readchomp(`git -C $repo rev-parse --short HEAD`)
    tables = ":!test/newcore/gates/measurements"
    clean = isempty(readchomp(`git -C $repo status --porcelain -- . $tables`))
    clean ? "commit $(head)" : "commit $(head) with uncommitted changes"
end
