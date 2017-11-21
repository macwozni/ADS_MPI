!> FRUIT unit test framework driver.

program fruit_test_driver
    use fruit
    use knot_vector_test

    call init_fruit()
    call init_fruit_xml()

    call test_knot1()
    call test_knot2()
    call test_knot3()

    call fruit_summary()
    call fruit_summary_xml()
    call fruit_finalize()

end program fruit_test_driver
