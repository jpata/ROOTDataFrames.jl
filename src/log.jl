const MAXDEPTH = 100
backtrace_list(bt) = [ccall(:jl_lookup_code_address, Any, (Ptr{Void}, Int32), b, 0) for b in bt]

IGNORED = [:anonymous, :eval_user_input, :rec_backtrace, :eval_body, :eval,:include,:include_from_node1,:process_options,:_start, :true_main]

macro debugprint(msg)
    return :(if DEBUG
        bt = backtrace()
        btl = backtrace_list(bt)
        names = Any[]
        for b in btl
            if (
                (b[1] in IGNORED) ||
                startswith(string(b[1]), "jl_") ||
                startswith(string(b[1]), "jl") ||
                startswith(string(b[1]), "julia")
            )
                continue
            end
            push!(names, string(b[1]))
        end
        println("[", join(reverse(names), ":"), "] " , $msg)
    end)
end
