is_delayed_array <- function(X) {
    methods::is(X, "DelayedMatrix") |
        methods::is(X, "DelayedArray")
}
