use theta_chars::points_in_ellipse;

#[test]
fn test_points_in_ellipse() {
    let r = 10000;
    assert_eq!(points_in_ellipse(0, 0, r).len(), 15701);
    assert_eq!(points_in_ellipse(0, 1, r).len(), 15708);
    assert_eq!(points_in_ellipse(1, 0, r).len(), 15684);
    assert_eq!(points_in_ellipse(1, 1, r).len(), 15716);
}
