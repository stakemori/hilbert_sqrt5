use theta_chars::points_in_ellipse;

#[test]
fn test_points_in_ellipse() {
    let r = 300;
    assert_eq!(points_in_ellipse(0, 0, r).len(), 63205);
    assert_eq!(points_in_ellipse(0, 1, r).len(), 63230);
    assert_eq!(points_in_ellipse(1, 0, r).len(), 63204);
    assert_eq!(points_in_ellipse(1, 1, r).len(), 63240);

    let r = 500;
    assert_eq!(points_in_ellipse(0, 0, r).len(), 175635);
    assert_eq!(points_in_ellipse(0, 1, r).len(), 175636);
    assert_eq!(points_in_ellipse(1, 0, r).len(), 175624);
    assert_eq!(points_in_ellipse(1, 1, r).len(), 175612);
}
