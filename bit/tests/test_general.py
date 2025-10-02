from unittest import mock
from bit.modules.general import tee

def test_tee_removes_ansi_and_logs(tmp_path):
    log_path = tmp_path / "test.log"
    msg_with_color = "\033[1mBold text\033[0m and normal"

    with mock.patch("sys.stdout", new_callable=mock.MagicMock()) as mock_stdout:
        tee(msg_with_color, log_path)

    # stdout write calls come as two separate writes: msg, then end
    calls = [mock.call(msg_with_color), mock.call("\n")]
    mock_stdout.write.assert_has_calls(calls)

    logged = log_path.read_text()
    assert "Bold text and normal" in logged
    assert "\033" not in logged

    log_path.unlink()
    tee("No newline", log_path, end="")
    assert log_path.read_text() == "No newline"
